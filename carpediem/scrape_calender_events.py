import argparse
import logging
import pandas as pd

from carpediem.calendar_utils import scrape_calendars

logging.basicConfig(level=logging.INFO)


def main(email_filename, calendar_directory, start_email, start_date, end_date, reprocess, limit):

    # Load list of emails from a csv
    df = pd.read_csv(email_filename)

    # Expects emails to be in the preferred_email or email column
    if "preferred_email" in df:
        emails = sorted(df["preferred_email"].combine_first(df["email"]).drop_duplicates())
    else:
        emails = sorted(df["email"].drop_duplicates())

    # Start from a given email if provided.
    if start_email is not None:
        start_pos = emails.index(start_email)
        emails = emails[start_pos:]

    # Scrape calendars
    scrape_calendars(
        emails,
        calendar_directory,
        reprocess=reprocess,
        start_date=start_date,
        end_date=end_date,
        limit=limit,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("email_filename", type=str)
    parser.add_argument("calendar_directory", type=str)
    parser.add_argument("--start_email", type=str, default=None)
    parser.add_argument("--start_date", type=str, default=None)
    parser.add_argument("--end_date", type=str, default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--reprocess", type=bool, default=True)
    args = parser.parse_args()

    main(
        args.email_filename,
        args.calendar_directory,
        args.start_email,
        args.start_date,
        args.end_date,
        args.reprocess,
        args.limit,
    )
