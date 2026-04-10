from __future__ import print_function
from dateutil.parser import parse
import json
import logging
import os.path
import pickle
import re

from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request

# If modifying these scopes, delete the file token.pickle.
SCOPES = ['https://www.googleapis.com/auth/calendar.readonly']


# Regex to identify rooms
# zymergen.com_2d3834373636353130313437@resource.calendar.google.com
ROOM_RE = re.compile('@resource.calendar.google.com')

# Regex to identify groups
# zymergen.com_2d3834373636353130313437@resource.calendar.google.com
GROUP_RE = re.compile('@group.calendar.google.com')

# Regex to identify real users
USER_RE = re.compile('[a-z]+@zymergen.com')


def create_service():
    creds = None
    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists('token.pickle'):
        with open('token.pickle', 'rb') as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file('credentials.json', SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open('token.pickle', 'wb') as token:
            pickle.dump(creds, token)

    service = build('calendar', 'v3', credentials=creds)
    return service


def list_calendars(service):
    calanders = []
    page_token = None

    while True:
        calendar_list = service.calendarList().list(
            pageToken=page_token, showHidden=True
        ).execute()
        calanders.extend(calendar_list.get('items', []))
        page_token = calendar_list.get('nextPageToken')
        if not page_token:
            break

    return calanders


def list_events(service, calendar_id, limit=None, start_date=None, end_date=None):
    events = []
    page_token = None

    if start_date is not None:
        time_min = parse(start_date).isoformat() + 'Z'  # 'Z' indicates UTC time
    else:
        time_min = None

    if end_date is not None:
        time_max = parse(end_date).isoformat() + 'Z'  # 'Z' indicates UTC time
    else:
        time_max = None

    while True:
        event_list = service.events().list(
            calendarId=calendar_id,
            singleEvents=True,
            pageToken=page_token,
            maxResults=limit or 1000,
            timeMin=time_min,
            timeMax=time_max,
            orderBy="startTime"
        ).execute()
        events.extend(event_list.get('items', []))
        page_token = event_list.get('nextPageToken')
        if not page_token:
            break
        if limit is not None and len(events) >= limit:
            events = events[:limit]
            break

    return events


def calendar_id_is_valid(service, calendar_id):
    try:
        service.events().list(
            calendarId=calendar_id, singleEvents=True, maxResults=1
        ).execute()
        return True
    except HttpError:
        return False


def is_room_id(rid):
    return bool(ROOM_RE.search(rid))


def is_group_id(rid):
    return bool(GROUP_RE.search(rid))


def get_event_room_ids(event):
    room_ids = set()
    for a in event.get('attendees', []):
        email = a["email"]
        if is_room_id(email):
            room_ids.add(email)
    return room_ids


def get_event_attendees(event):
    attendees = [
        a for a in event.get('attendees', [])
        if a.get('responseStatus') != "declined"
        and not is_room_id(a["email"])
        and not is_group_id(a["email"])
        and USER_RE.search(a["email"])
    ]

    return [a["email"].split("@")[0] for a in attendees]


def scrape_calendars(calendar_ids, calendar_directory, reprocess=None, **kwargs):
    # Create calendar service
    service = create_service()

    # Iterate over calendars.
    for cid in calendar_ids:
        cname = cid.split("@")[0]
        events_filename = os.path.join(calendar_directory, cname + ".json")

        # Optionally skip calendars that already have data.
        if not reprocess and os.path.isfile(events_filename):
            logging.info(f"Skipping calendar event scraping for calendar: {cname}")
            continue

        # Check if the user id is valid.
        logging.info(f"Checking for valid calendar: {cname}")
        if calendar_id_is_valid(service, cid):

            # Scrape calendar events
            logging.info(f"Scraping events for calendar: {cname}")
            events = list_events(service, cid, **kwargs)

            # Save results as json
            logging.info(f"Saving calendar events to {events_filename}")
            with open(events_filename, "w") as f:
                json.dump(events, f)

        else:
            logging.warning(f"Unable to locate calendar: {cname}.")


def iter_calendars(calendar_directory):
    directory = os.fsencode(calendar_directory)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".json"):
            with open(os.path.join(calendar_directory, filename), 'r') as f:
                yield json.load(f)


def iter_events(events, min_attendees=2, max_attendees=10, skip_repeat_events=True):
    for event in events:
        attendees = get_event_attendees(event)
        if not (min_attendees <= len(attendees) <= max_attendees):
            continue
        if skip_repeat_events and event.get('recurringEventId'):
            continue
        yield event


def get_room_ids(calendar_directory):
    room_ids = set()
    for events in iter_calendars(calendar_directory):
        for event in iter_events(events):
            room_ids.update(get_event_room_ids(event))
    return sorted(list(room_ids))
