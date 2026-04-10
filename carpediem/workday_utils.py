import logging
import re


def full_name_to_email(full_name):
    # Remove non alpha numeric
    full_name = re.sub('\W+ ', '', full_name)
    # Remove I, II, etc
    full_name = re.sub(' II+', '', full_name)
    # Remove -_
    full_name = re.sub('[-_]', ' ', full_name)
    # Remove weird contract label
    full_name = full_name.replace(" [C]", "")
    # Split name into parts
    split_name = full_name.lower().split(' ')
    if len(split_name) > 2:
        logging.warning(
            f"Full name contains more than two names: {full_name}. "
            "The generated email may not be correct."
        )
    return split_name[0][0] + ''.join(split_name[1:]) + "@zymergen.com"
