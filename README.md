# Ligify
This repo focuses on predicting protein-ligand interactions.

# Purpose

Following the successful release of a web tool for Ligify, the development team moved into creating a publicly available service which is what this code contains. 

# Running

In order to run Ligify locally, you need the following:

- [aws-sam-cli](https://github.com/aws/aws-sam-cli)
- [Visual Studio Code](https://code.visualstudio.com/download) or your favorite text editor

# Environment

In the root of the project, create a `env.local.json` file with the following structure:

```
{
    "Parameters": {
        "NCBI_API_KEY": "XXXXXXXXXXXXXXX"
    }
}
```

And NCBI API Key can be obtained by registering for an [NCBI Account](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us) and requesting an API key. Without this, you'll need to modify the code slightly and understand the [rate limits](https://support.nlm.nih.gov/knowledgebase/article/KA-05318/en-us) associated without a key. 

# Running

In order to run Ligify locally, you have two options:

1. Open two terminals side by side and run the following commands:

- `sam build --use-container` <- terminal 1
- `sam local start-api --env-vars env.local.json` <- terminal 2 after step #1 completes

2. Add [nodemon](https://www.npmjs.com/package/nodemon) globally, [concurrently](https://www.npmjs.com/package/concurrently) globally and python3.11 then run the following command:

- `sam build && concurrently "nodemon --on-change-only --ext py --exec sam build" "sam local start-api --env-vars env.local.json"`

Refer to issue [#921](https://github.com/aws/aws-sam-cli/issues/921) for more information

# Compromises

During development, the team came to realize that numpy did not behave properly in AWS SAM when using just `local start-api`. See [#7429](https://github.com/aws/aws-sam-cli/issues/7429) and [#27338](https://github.com/numpy/numpy/issues/27338) for more information. To get around this, the `--use-container` flag was introduced which creates an isolated environment locally and resolves dependency issues.

Occasionally, the rdkit wheel fails to build. When this happens, simply run the your appropriate build command again (mostly seems to happen with concurrently)