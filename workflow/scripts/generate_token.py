import os
import sys
import yaml
import argparse
import getpass
import requests


def main():
    parser = argparse.ArgumentParser()
    argument_parser(parser)
    args = parser.parse_args()
    with open(args.config, "r") as config_file:
        config = yaml.safe_load(config_file)
    user = input("Enter SHIP user:")
    password = getpass.getpass("Enter password:")
    try:
        request = requests.get(
            config["ship"]["generate_token_url"], auth=(user, password)
        )
    except:
        print(
            "Requesting token failed. Please validate url, user and password.",
            file=sys.stderr,
        )
        sys.exit(1)
    if request.ok:
        with open(config["ship"]["token_path"], "w") as f:
            print(request.text.strip(), file=f)
    else:
        raise ValueError(
            f"Requesting token failed returning: {request.status_code} - {request.text}"
        )


def argument_parser(parser):
    parser.add_argument(
        "-c",
        "--config",
        help="Path of configuration file as yaml",
        type=str,
        required=True,
    )


if __name__ == "__main__":
    main()
