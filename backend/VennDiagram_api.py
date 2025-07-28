import requests

# Define the endpoint and parameters
url = "http://localhost:8000/hello"
params = {"name": "Ryan"}

# Send GET request
response = requests.get(url, params=params)

# Show the response
print(response.json())
