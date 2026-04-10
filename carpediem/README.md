# 📅 Carpe Diem

> *Seize the day — and all of its calendar data.*

A 2019 hackday project that mines Google Calendar events across an organization to build a
**social graph of who meets with whom**. The output is a weighted network you can visualize in
[Cytoscape.js](https://js.cytoscape.org/) (or any other graph tool that accepts the format).
Graph analysis offers perspective on organizational design and individual influence.

---

## How it works

The pipeline has three stages:

1. **Process the org chart** — takes a CSV export from Workday (or similar HR tools) and converts
   employee display names into email addresses.
2. **Scrape calendars** — uses the Google Calendar API to pull meeting events for every employee and
   saves them locally as JSON.
3. **Build the graph** — iterates over all scraped events, and for each meeting adds (or increments)
   a weighted edge between every pair of co-attendees. The result is exported as a NetworkX pickle
   and a Cytoscape.js JSON file.

---

## Project structure

```
carpediem/
├── __init__.py
├── calendar_utils.py          # Google Calendar API helpers, graph building logic
├── workday_utils.py           # Name → email conversion
├── process_org_chart.py       # CLI: enrich an org chart CSV with email addresses
├── scrape_calender_events.py  # CLI: scrape calendar events for a list of employees
├── process_calendar_events.py # CLI: build the co-meeting graph from scraped data
└── demo.py                    # Quick sanity check — prints your next 10 calendar events
```

---

## Setup

### Prerequisites

- Python 3.6+
- A Google Cloud project with the **Google Calendar API** enabled
- A `credentials.json` file from the Google Cloud Console (OAuth 2.0 client)

### Install dependencies

```bash
pip install google-api-python-client google-auth-oauthlib google-auth-httplib2 \
            python-dateutil networkx pandas
```

### Authenticate

On first run, the script will open a browser window to complete OAuth. A `token.pickle` file will be
saved locally so you won't need to log in again.

---

## Usage

### 1. Process the org chart

Convert a Workday/HR CSV export (with a `Name` column) into email addresses:

```bash
python -m carpediem.process_org_chart org_chart.csv
```

This adds an `email` column to the CSV in-place.

### 2. Scrape calendar events

Pull meeting events for every employee in the CSV:

```bash
python -m carpediem.scrape_calender_events org_chart.csv ./calendar_data/ \
    --start_date 2019-01-01 \
    --end_date   2019-09-01
```

Events are saved as one JSON file per person inside `./calendar_data/`. Already-scraped calendars
are skipped by default.

### 3. Build the co-meeting graph

```bash
# Build the weighted social graph
python -m carpediem.process_calendar_events ./calendar_data/ build_graph

# (Optional) extract all room resource IDs found in the data
python -m carpediem.process_calendar_events ./calendar_data/ room_ids
```

**Outputs:**
- `network.gpickle` — NetworkX graph for further analysis in Python
- `network.cyjs` — Cytoscape.js JSON for interactive visualization

---

## Graph details

- **Nodes** — individual employees (identified by username)
- **Edges** — weighted by co-meeting frequency (each shared meeting adds `0.5` to the edge weight)
- **Filtered out:** room resources, Google Group calendars, declined invites, recurring events, and
  meetings with fewer than 2 or more than 10 attendees

---


## What you can do with the graph

The weighted edges let you run standard organizational network analysis: find **central connectors**
(high betweenness centrality) who bridge otherwise-disconnected parts of the org, spot **natural
clusters** that don't match the official org chart, identify **cross-team weak ties**, or flag
**isolates** who barely appear. Edge weight distinguishes a shared all-hands from a close working
relationship.

---

## Gotchas

### Private calendars cause asymmetric undercounting
No admin access is needed — a standard org-member OAuth token can read any calendar shared org-wide,
which is the default in most Google Workspace orgs. People with private calendars aren't completely
invisible either: since `get_event_attendees` reads the full attendee list from each event, a
private-calendar person still shows up via meetings scraped from their colleagues' public
calendars. What you lose is any meeting they hosted or attended where every other participant also
had a private calendar. The net effect is that private-calendar people tend to look slightly less
connected than they are — treat unusually sparse nodes with some skepticism.

### The graph reflects one-off meetings, not routines
`iter_events` skips recurring events by design. Weekly 1:1s and standups would otherwise dominate
and drown out ad-hoc collaboration signals. Worth keeping in mind when interpreting results.

### Each event appears in multiple scraped files
Since calendars are scraped per-person, a shared meeting lands in every attendee's JSON. The
`processed_events` set in `build_graph` deduplicates by event ID — don't remove it.

### `--reprocess` defaults to `True`, defeating the cache
The CLI will re-scrape every calendar from scratch unless you pass `--reprocess False`. This is
almost certainly a bug — the caching logic exists precisely so interrupted scrapes can resume.

### Rate limiting
No retry or backoff logic. For large orgs, add a `time.sleep()` between requests or handle
`HttpError 429`.

### `nx.write_gpickle` was removed in NetworkX 3.0
Replace it with:
```python
import pickle
with open("network.gpickle", "wb") as f:
    pickle.dump(G, f)
```

---


## Notes

- The Zymergen-specific email domain (`@zymergen.com`) and regex patterns are hard-coded in
  `calendar_utils.py` — swap these out for your own org's domain.
- This was a one-day hackday project. The code works but is not production-hardened — error handling
  is minimal and there are a few typos in filenames (`scrape_calender_events.py`) that have been
  preserved for posterity.
- `list_calendars()` in `calendar_utils.py` is defined but never called by the pipeline — it was likely useful for exploration during the hackday.
