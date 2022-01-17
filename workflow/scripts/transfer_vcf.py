import requests

# https://stackoverflow.com/questions/43246857/python-requests-post-a-file
# https://stackoverflow.com/questions/68477/send-file-using-post-from-a-python-script

token_path = "config/ship_token.txt"
endpoint_url = (
    ""
)
vcf = ""
vcf_name = ""

with open(token_path, "r+") as f:
    token = f.read().strip()
headers = {"Authorization": "Bearer " + token}
print(headers)
with open(vcf, "rb") as vcf_file:
    files = {"filename": vcf_name, "content": vcf_file}
    response = requests.post(endpoint_url, headers=headers, files=files)
if response.ok:
    print(response.json())
else:
    raise ValueError(f"Patient lookup failed: {response.status_code} - {response.text}")
