import requests


def update_token():
    with open(config["ship"]["token_path"], "r+") as f:
        token = f.read().strip()
        headers = {"Authorization": "Bearer " + token}
        print(headers)
        response = requests.get(config["ship"]["refresh_token_url"], headers=headers)
        if response.ok:
            new_token = response.text
        else:
            raise ValueError(
                f"Updating token failed: {response.status_code} - {response.text}"
            )
        f.seek(0)
        f.write(new_token)
        f.truncate()
