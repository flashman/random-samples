import argparse
import pandas as pd

from carpediem.workday_utils import full_name_to_email


def main(filename):
    df = pd.read_csv(filename)
    df["email"] = df['Name'].apply(full_name_to_email)
    df.to_csv(filename, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()

    main(args.filename)
