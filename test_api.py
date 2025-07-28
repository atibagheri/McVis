import requests
import json

url = "http://localhost:8000/venn-upset"
headers = {"Content-Type": "application/json"}
data = {"paths": ["../Ex_expression1.txt", 
        "../Ex_expression2.txt", 
        "../Ex_expression6.txt", 
        "../Ex_expression7.txt"]
        }

# try:
response = requests.post(url, headers=headers, data=json.dumps(data))
print(response.content.output_path)
