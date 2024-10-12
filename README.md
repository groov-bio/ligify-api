# Ligify
This repo focuses on predicting protein-ligand interactions and providing a publically available service to work with this prediction algorithm.

# Purpose

Following the successful release of a web tool for Ligify, the development team moved into creating a publicly available service which is what this code contains. 

# Running

In order to run Ligify locally, you need the following:

- [docker](https://docs.docker.com/engine/install/)
- [aws-sam-cli](https://github.com/aws/aws-sam-cli)
- [Visual Studio Code](https://code.visualstudio.com/download) or your favorite text editor

# Environment

In the root of the project, create a `env.local.json` file with the following structure:

```
{
    "Parameters": {
        "NcbiApiKey": "XXXXXXXXXXXXXXX"
    }
}
```

And NCBI API Key can be obtained by registering for an [NCBI Account](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us) and requesting an API key. Without this, you'll need to modify the code slightly and understand the [rate limits](https://support.nlm.nih.gov/knowledgebase/article/KA-05318/en-us) associated without a key. 

# Running

In order to run Ligify locally, you have two options:

### Option A

Open two terminals side by side and run the following commands:

- `sam build -t local.template.yaml` <- `terminal #1`
- `sam local start-api --env-vars env.local.json` <- `terminal #2 after step #1 completes`

### Option B

Add [nodemon](https://www.npmjs.com/package/nodemon) globally, [concurrently](https://www.npmjs.com/package/concurrently) globally and python3.11 then run the following command:

- `sam build -t local.template.yaml && concurrently "nodemon --on-change-only --ext py --exec sam build -t local.template.yaml" "sam local start-api --env-vars env.local.json"`

Refer to issue [#921](https://github.com/aws/aws-sam-cli/issues/921) for more information

You'll know it's ready to use locally when the output contains `Running on http://127.0.0.1:3000`

# Deployment

The following command can be used to deploy your own cloudformation stack:

`sam deploy --stack-name ligify-api --resolve-s3 --capabilities CAPABILITY_IAM --resolve-image-repos --guided --debug`

# Venv

It's recommended to install the packages in `/ligify` locally using a virtual environment. Depending on your system, installing can be different but the [venv](https://docs.python.org/3/library/venv.html) docs generally cover most systems. The directory `ligify-venv` is already part of `.gitignore` so it's suggested to use that naming convention.

# Contributing

We welcome any pull requests to help make the ligify-api better. Review the below code formatting guidelines and be sure to fill in the PR template appropriately.

# Code formatting and quality

We use [ruff](https://github.com/astral-sh/ruff) for python code formatting and quality which is installed independently of the bundled dependencies. Please review their documentation in order to properly format your codebase. 

You can run `ruff check` for quality and `ruff format` for formatting

# Compromises

1. During development, the team came to realize that numpy did not behave properly in AWS SAM when using just `local start-api`. See [#7429](https://github.com/aws/aws-sam-cli/issues/7429) and [#27338](https://github.com/numpy/numpy/issues/27338) for more information. To get around this, the `--use-container` flag was introduced which creates an isolated environment locally and resolves dependency issues.

2. Occasionally, the rdkit wheel fails to build. When this happens, simply run the your appropriate build command again (mostly seems to happen with concurrently).

3. Cloudfront has a hard limit of 180 seconds which you need to request a quota increase for. In the event of a 504 timeout, you should default the request back to the function URL directly.

4. In order to add a WebACL, you need to deploy your stack into us-east-1. 