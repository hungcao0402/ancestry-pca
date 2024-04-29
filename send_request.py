import requests
import json
import base64
import io
import matplotlib.pyplot as plt
from PIL import Image

def send_prediction_request(population, superpopulation, user):
    url = 'http://0.0.0.0:8080/v1/models/AncestryPCA:predict'  # Replace with your actual endpoint URL
    
    # Prepare the request payload
    payload = {
        "instances": [
            {
                "population": population,
                "superpopulation": superpopulation,
                "user": user
            }
        ]
    }
    
    # Send POST request with JSON payload
    response = requests.post(url, json=payload)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Parse the response JSON
        predictions = response.json()
        return predictions
    else:
        print(f"Request failed with status code: {response.status_code}")
        return None

# Example usage
population = None
superpopulation = "EUR EAS AFR SAS AMR"
user = "HG01595 HG01607 HG03702 HG03078"
predictions = send_prediction_request(population, superpopulation, user)



image_data = base64.b64decode(predictions['predictions']['b64'])

# Convert the decoded image data into a BytesIO object
image_buffer = io.BytesIO(image_data)

# Read the image from the BytesIO object using matplotlib
# Save the figure to a file
im = Image.open(image_buffer)
im.save('temp.png', 'PNG')

