FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.138.cpg1-1

ENV PYTHONDONTWRITEBYTECODE=1
ENV VERSION=0.1.1

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /workflow_name

COPY src src/
COPY LICENSE pyproject.toml README.md ./

RUN pip install --no-cache-dir .
