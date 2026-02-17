# Somatic Variant Classifier --- Task II

## Containerization & Workflow Orchestration with k3d + Argo Workflows

------------------------------------------------------------------------

## Overview

This project containerizes a somatic variant classifier (Task I) and
deploys it as a reusable workflow on a multi-node local Kubernetes
cluster using:

-   Docker (containerization)
-   k3d (local multi-node Kubernetes cluster)
-   Argo Workflows (workflow orchestration)

The system allows parametrized execution of the classifier on tumor and
control VCF samples in a reproducible and portable manner.

------------------------------------------------------------------------

# Project Structure

Application source code and large reference datasets are intentionally
separated.

## 📁 Application Repository (version controlled)

    ~/Desktop/somatic-classifier-project/
    ├── Dockerfile
    ├── requirements.txt
    ├── workflowtemplate.yaml
    └── germline_somatic_classifier/
        └── somatic_variant_classifier.py

## 📁 External Reference Data (NOT committed to Git)

    ~/refdata/
    ├── gnomad4_1pct.vcf.gz
    ├── gnomad4_1pct.vcf.gz.tbi
    ├── hg38_cosmic91.txt.gz
    ├── hg38_cosmic91.txt.gz.tbi
    ├── samples/
    └── results/

Large genomic reference datasets are intentionally excluded from Git for
ethical, size, and reproducibility reasons.

------------------------------------------------------------------------

# Reference File Preparation (IMPORTANT)

The classifier requires indexed reference files (`.tbi`).

If index files are missing, the workflow will fail.

Install tabix:

``` bash
sudo apt-get update
sudo apt-get install -y tabix
```

### gnomAD indexing

``` bash
tabix -p vcf ~/refdata/gnomad4_1pct.vcf.gz
```

### COSMIC indexing

If starting from plain `.txt`:

``` bash
bgzip ~/refdata/hg38_cosmic91.txt
tabix -s 1 -b 2 -e 2 ~/refdata/hg38_cosmic91.txt.gz
```

If already gzipped but not indexed:

``` bash
tabix -s 1 -b 2 -e 2 ~/refdata/hg38_cosmic91.txt.gz
```

Indexed access enables streaming queries and prevents OOM errors,
especially in low‑memory environments.

------------------------------------------------------------------------

# Architecture

    Local Machine
    └── k3d Kubernetes Cluster
        ├── 1 Server Node
        ├── 3 Agent Nodes
        └── Argo Workflows
            └── WorkflowTemplate
                └── Runs classifier container

Reference data is mounted from the host machine into all cluster nodes:

``` bash
k3d cluster create somatic-cluster \
  --agents 3 \
  --volume ~/refdata:/refdata@all
```

------------------------------------------------------------------------

# Tested Environment

Successfully validated on:

-   OS: Ubuntu 24.04 (Linux kernel 6.14.x)
-   Architecture: x86_64
-   RAM: 4 GB
-   Docker 20+
-   kubectl 1.30+
-   k3d 5.x
-   Argo CLI 3.5.x

The limited 4GB RAM environment motivated memory‑efficient streaming and
indexed file access.

------------------------------------------------------------------------

# Step-by-Step Setup

## 1️⃣ Create k3d Cluster

``` bash
k3d cluster create somatic-cluster --agents 3 --volume ~/refdata:/refdata@all
kubectl get nodes
```

### Expected Output:

``` bash
NAME                           STATUS   ROLES                  AGE   VERSION
k3d-somatic-cluster-agent-0    Ready    <none>                 23s   v1.31.5+k3s1
k3d-somatic-cluster-agent-1    Ready    <none>                 25s   v1.31.5+k3s1
k3d-somatic-cluster-agent-2    Ready    <none>                 24s   v1.31.5+k3s1
k3d-somatic-cluster-server-0   Ready    control-plane,master   55s   v1.31.5+k3s1
```


## 2️⃣ Install Argo Workflows (v3.5.11)

``` bash
kubectl create namespace argo

kubectl apply -n argo \
  -f https://github.com/argoproj/argo-workflows/releases/download/v3.5.11/install.yaml
```

Verification:

``` bash
kubectl get pods -n argo -w
```

Expected Output:

``` bash
NAME                                   READY   STATUS    RESTARTS   AGE
argo-server-59df945c7-5b6fq            1/1     Running   0          4m10s
workflow-controller-788b597c74-tknn2   1/1     Running   0          4m9s

```

Patch for local development:

``` bash
kubectl patch deployment argo-server -n argo \
  --type='json' \
  -p='[{"op":"replace","path":"/spec/template/spec/containers/0/args","value":["server","--auth-mode=server"]}]'

kubectl rollout restart deployment argo-server -n argo
```

Access UI:

``` bash
kubectl port-forward svc/argo-server -n argo 2746:2746
```

Open:

    https://localhost:2746

------------------------------------------------------------------------

## 3️⃣ Build Docker Image

``` bash
docker build -t somatic-classifier:latest .
k3d image import somatic-classifier:latest -c somatic-cluster
```

Standalone verification:

``` bash
docker run --rm \
  -v ~/refdata:/refdata \
  somatic-classifier:latest \
  --sample-vcf /refdata/samples/CO8-PA-26_mutect2_vlod.vcf.gz \
  --gnomad-vcf /refdata/gnomad4_1pct.vcf.gz \
  --cosmic-tsv /refdata/hg38_cosmic91.txt.gz \
  --output /refdata/results/test_output.tsv
```

Exit code must be 0 on success.

Expected Output:
``` bash
Starting Somatic vs Germline Variant Classifier
Completed successfully in 5.58 seconds. Exit Code: 0
```
------------------------------------------------------------------------

## 4️⃣ Apply WorkflowTemplate

``` bash
kubectl apply -f workflowtemplate.yaml
```

Template properties:

-   Accepts parameter: `sample-vcf`
-   Mounts `/refdata`
-   Writes results to `/refdata/results`
-   Uses `serviceAccountName: argo`

------------------------------------------------------------------------

## 5️⃣ Submit Workflow

Tumor sample:

``` bash
argo submit -n argo \
  --from workflowtemplate/somatic-classifier-template \
  -p sample-vcf=CO8-PA-26_mutect2_vlod.vcf.gz \
  --watch
```

Example Output:
``` bash
Name:                somatic-classifier-template-j49jr
Namespace:           argo
ServiceAccount:      argo
Status:              Succeeded
Conditions:          
 PodRunning          False
 Completed           True
Created:             Tue Feb 17 23:46:29 +0300 (46 seconds ago)
Started:             Tue Feb 17 23:46:29 +0300 (46 seconds ago)
Finished:            Tue Feb 17 23:47:15 +0300 (now)
Duration:            46 seconds
Progress:            1/1
ResourcesDuration:   2s*(1 cpu),31s*(100Mi memory)
Parameters:          
  sample-vcf:        CO8-PA-26_mutect2_vlod.vcf.gz

STEP                                  TEMPLATE        PODNAME                            DURATION  MESSAGE
 ✔ somatic-classifier-template-j49jr  run-classifier  somatic-classifier-template-j49jr  36s 
```

Control sample:

``` bash
argo submit -n argo \
  --from workflowtemplate/somatic-classifier-template \
  -p sample-vcf=CO8-NC-pop_mutect2_vlod.vcf.gz \
  --watch
```

Monitoring:

``` bash
argo list -n argo
argo get <workflow-name> -n argo
argo logs <workflow-name> -n argo
```

------------------------------------------------------------------------

# Output Verification

Results are written to:

    ~/refdata/results/

Example:

    somatic-classifier-template-xxxxx.tsv

Each workflow execution generates deterministic TSV output.

------------------------------------------------------------------------

# Known Issues & Resolutions

✔ Argo v4 CRD size error → Downgraded to v3.5.11\
✔ TLS readiness probe issue → Patched auth mode\
✔ ServiceAccount permission error → Added `serviceAccountName: argo`\
✔ Artifact storage error → Removed artifact configuration and used
host-mounted storage

All issues resolved.

------------------------------------------------------------------------

# Task II Compliance Summary

✔ Multi-node k3d cluster (1 server + 3 agents)\
✔ Argo Workflows operational\
✔ Dockerized classifier\
✔ Standalone container validation\
✔ Parameterized WorkflowTemplate\
✔ Indexed reference datasets\
✔ Tumor & control workflow executions successful\
✔ Deterministic outputs generated\
✔ No external artifact repository required
