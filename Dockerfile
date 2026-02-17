FROM python:3.10-slim

WORKDIR /app

# System deps for cyvcf2 / pysam
RUN apt-get update && apt-get install -y \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    tabix \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY germline_somatic_classifier/ ./germline_somatic_classifier/

ENTRYPOINT ["python", "-m", "germline_somatic_classifier.somatic_variant_classifier"]
