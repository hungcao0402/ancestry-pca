# Use the official lightweight Python image.
FROM python:3.11-slim
 
ENV APP_HOME /app
WORKDIR $APP_HOME

# Install production dependencies.
COPY requirements.txt ./

COPY ./ ./

RUN pip install --no-cache-dir -r ./requirements.txt
 
# Copy local code to the container image
# COPY kf_model_server.py ./
 
CMD ["python", "kf_model_server.py"]
