from typing import List, Dict
from kserve import Model, ModelServer
from ancestry_pca import ancestry_pca
import base64
import io 
import matplotlib.pyplot as plt

# Define the custom model server class
class AncestryPCA(Model):
    def __init__(self, name: str):
       super().__init__(name)
       self.name = name
       self.load()


    def load(self):
        self.ready = True
        

    def predict(self, request: Dict, headers: Dict[str, str] = None) -> Dict:
        population = request["instances"][0]['population']
        superpopulation = request["instances"][0]['superpopulation']
        user = request["instances"][0]['user']

        population = None if population == '' else population
        superpopulation = None if superpopulation == '' else superpopulation
        user = None if user == '' else user

        figure = ancestry_pca(p=population, sp=superpopulation, user=user)

        buffer = io.BytesIO()
        figure.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
        buffer.seek(0)

        # Encode the image data to Base64
        b64 = base64.b64encode(buffer.getvalue()).decode('utf-8')

        # Close the buffer
        buffer.close()

        return {"predictions": {'b64': b64}
        }

if __name__ == "__main__":
    model = AncestryPCA("AncestryPCA")
    ModelServer().start([model])
