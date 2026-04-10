import argparse
import json
import logging
import networkx as nx
import pandas as pd

from carpediem.calendar_utils import iter_calendars, iter_events, get_event_attendees, get_room_ids

logging.basicConfig(level=logging.INFO)


def extract_room_ids(calendar_directory):
    room_ids = get_room_ids(calendar_directory)
    df = pd.DataFrame({"email": room_ids})
    df.to_csv("room_ids.csv", index=False)


def build_graph(calendar_directory):
    # Construct multi-graph
    processed_events = set()
    G = nx.Graph()
    for events in iter_calendars(calendar_directory):
        for event in iter_events(events):
            eid = event['id']
            if eid not in processed_events:
                attendees = get_event_attendees(event)
                for a1 in attendees:
                    for a2 in attendees:
                        if a1 != a2:
                            if not G.has_edge(a1, a2):
                                G.add_edge(a1, a2, weight=0.5)
                            else:
                                G[a1][a2]['weight'] += 0.5
                processed_events.add(eid)
            else:
                logging.info(f'Event {eid} has already been processed')

    nx.write_gpickle(G, "network.gpickle")
    with open('network.cyjs', 'w') as f:
        json.dump(nx.cytoscape_data(G), f)


def main(calendar_directory, mode):
    if mode == "room_ids":
        extract_room_ids(calendar_directory)
    elif mode == "build_graph":
        build_graph(calendar_directory)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("calendar_directory", type=str)
    parser.add_argument("mode", type=str, choices=["room_ids", "build_graph"])
    args = parser.parse_args()

    main(
        args.calendar_directory, args.mode
    )
