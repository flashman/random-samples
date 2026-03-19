import os
import sys

test_dir = os.path.dirname(os.path.realpath(__file__))

os.environ["DB_USERNAME_PATH"] = os.path.join(test_dir, "db_username.txt")
os.environ["DB_PASSWORD_PATH"] = os.path.join(test_dir, "db_password.txt")


sys.path.append(os.path.abspath("../.code"))
