import re
import json
import time
import csv
import os
import threading
from collections import deque
from dotenv import load_dotenv

import boto3
from boto3.dynamodb.conditions import Key
from marshmallow import Schema, fields, ValidationError, validate

import requests

from fetch_data import fetch_data
from genbank.create_genbank import create_genbank
from predict.pubchem import get_inchikey

##########################################################
# VALIDATION SCHEMAS
##########################################################

def validate_bool(value):
    if not isinstance(value, bool):
        raise ValidationError("reviewed requires a boolean value.")


class FilterSchema(Schema):
    max_reactions = fields.Integer(
        required=True,
        validate=validate.Range(
            min=1, error="max_reaction must be a positive integer."
        ),
    )
    proteins_per_reaction = fields.Integer(
        required=True,
        validate=validate.Range(
            min=1, error="proteins_per_reaction must be a positive integer."
        ),
    )
    reviewed = fields.Boolean(required=True, validate=validate_bool)
    lineage = fields.Str(
        required=True,
        validate=validate.OneOf(
            ["Phylum", "Class", "Order", "Family", "Genus", "None"],
            error="lineage should be one of Phylum, Class, Order, Family, Genus, None",
        ),
    )
    max_operons = fields.Integer(
        required=True,
        validate=validate.Range(
            min=1, error="max_operons must be a positive integer."
        ),
    )
    max_alt_chems = fields.Integer(
        required=True,
        validate=validate.Range(
            min=1, error="max_alt_chems must be a positive integer."
        ),
    )


class InputSchema(Schema):
    # This might need to be adjusted in the future
    # https://gist.github.com/lsauer/1312860/63687490915641fdb4e377cf581211da9d5d64c0
    smiles = fields.String(
        validate=validate.Regexp(
            r"^([^J][a-z0-9@+\-\[\]\(\)\\\/%=#$]{6,})$",
            flags=re.IGNORECASE,
            error="Invalid SMILES provided",
        )
    )
    filters = fields.Nested(FilterSchema)

##########################################################
# DYNAMODB SETUP
##########################################################

def ensure_table_exists():
    try:
        # Check if table exists
        table = dynamodb.Table("Chemicals")
        table.table_status
        print("Table exists:", table.table_name)
        return table
    except Exception as e:
        print("Table doesn't exist, creating...")
        try:
            # Create the table
            table = dynamodb.create_table(
                TableName='Chemicals',
                KeySchema=[
                    {
                        'AttributeName': 'SMILES',
                        'KeyType': 'HASH'  # Partition key
                    },
                    {
                        'AttributeName': 'chunk_index',
                        'KeyType': 'RANGE'  # Sort key
                    }
                ],
                AttributeDefinitions=[
                    {
                        'AttributeName': 'SMILES',
                        'AttributeType': 'S'
                    },
                    {
                        'AttributeName': 'chunk_index',
                        'AttributeType': 'N'
                    }
                ],
                ProvisionedThroughput={
                    'ReadCapacityUnits': 5,
                    'WriteCapacityUnits': 5
                }
            )
            # Wait until the table exists
            table.meta.client.get_waiter('table_exists').wait(TableName='Chemicals')
            print("Table created successfully!")
            return table
        except Exception as create_error:
            print("Error creating table:", create_error)
            raise

# Connect to local DynamoDB
dynamodb = boto3.resource(
    'dynamodb',
    endpoint_url="http://host.docker.internal:8000",
    region_name="localhost", # region arbitrary since local
    aws_access_key_id="fakeMyKeyId",
    aws_secret_access_key="fakeSecretAccessKey"
)
table = ensure_table_exists()

def store_in_dynamodb(smiles, data, data_is_error = False):
    """
    Store the given data under the item keyed by SMILES.
    If data is >350KB once marshalled, break it into smaller pieces.
    We'll store them as chunked items with a sort key or suffix.
    """
    # Convert data to JSON string
    data_str = json.dumps(data)
    data_size = len(data_str.encode('utf-8'))  # approximate size in bytes

    max_chunk_size = 350000  # 350 KB approx
    if data_size <= max_chunk_size:
        # Just one write
        table.put_item(
            Item={
                "SMILES": smiles,
                "chunk_index": 0,
                "data": data if not data_is_error else None,
                "error": data if data_is_error else None
            }
        )
    else:
        # Break into chunks
        chunks = []
        start = 0
        index = 0
        while start < data_size:
            end = start + max_chunk_size
            chunk_data = data_str[start:end]
            chunks.append(chunk_data)
            start = end
        # Put each chunk with a separate chunk_index
        for i, chunk in enumerate(chunks):
            table.put_item(
                Item={
                    "SMILES": smiles,
                    "chunk_index": i,
                    "data": chunk
                }
            )

##########################################################
# MAIN PROCESSING LOGIC
##########################################################

def create_plasmid(regulators, chemical):
    for regulator in regulators:
        result = create_genbank(
            regulator["refseq"],
            chemical,
            regulator["protein"]["context"]["promoter"]["regulated_seq"],
            regulator["reg_protein_seq"],
        )
        regulator["plasmid_sequence"] = str(result)
    return regulators

def process_batch(filters):
    # Read merged_db.csv, take first 5 entries
    csv_path = "./merge_db.csv"
    if not os.path.exists(csv_path):
        print("merge_db.csv not found!")
        return

    results = []
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        entries = list(reader)[:100]  # first 5 entries
        for entry in entries:
            smiles = entry["SMILES"]
            chemical_name = entry["Name"]
            chebid = entry["CHEBID"]

            # We'll construct the full result only if successful
            # On error, we store only SMILES and error message
            try:
                # Attempt to get the InChIKey
                try:
                    inchi_key = get_inchikey(smiles, "smiles")
                except Exception as e:
                    error_item = {
                        "error": f"Failed to get InChIKey: {str(e)}"
                    }
                    store_in_dynamodb(smiles, error_item, True)
                    results.append(error_item)
                    continue  # Move on to the next entry

                # Attempt to fetch data
                try:
                    regulators, metrics = fetch_data(inchi_key, filters)
                except Exception as e:
                    error_item = {
                        "error": f"Failed to fetch data: {str(e)}"
                    }
                    store_in_dynamodb(smiles, error_item, True)
                    results.append(error_item)
                    continue

                # Attempt to create plasmids
                try:
                    regulators = create_plasmid(regulators, chemical_name)
                except Exception as e:
                    error_item = {
                        "error": f"Failed to create plasmid: {str(e)}"
                    }
                    store_in_dynamodb(smiles, error_item, True)
                    results.append(error_item)
                    continue

                # If everything succeeded, store the full data
                full_item = {
                    "CHEBID": chebid,
                    "SMILES": smiles,
                    "Name": chemical_name,
                    "metrics": metrics,
                    "regulators": regulators
                }
                store_in_dynamodb(smiles, full_item)
                results.append(full_item)

            except Exception as e:
                # Catch any other unforeseen exceptions at the outer level
                error_item = {
                    "SMILES": smiles,
                    "error": f"Unexpected error: {str(e)}"
                }
                store_in_dynamodb(smiles, error_item)
                results.append(error_item)
                # Continue with the next entry

    # After batch process completes, just return
    return results


##########################################################
# LAMBDA HANDLER
##########################################################

def lambda_handler(event, context):
    path = event.get('rawPath') or event.get('path')
    method = event.get('httpMethod', '').upper()
    print(f"Method: {method}, Path: {path}")  # Log the method and path

    if not path or '/ligify' not in path:
        return generate_response(403, "Forbidden")

    if method == 'OPTIONS':
        return generate_response(200, "", is_options=True)

    load_dotenv()

    try:
        body = json.loads(event.get("body", "{}"))
    except json.JSONDecodeError:
        return generate_response(
            400,
            {"message": "Invalid JSON in request body."}
        )

    input_schema = InputSchema()

    try:
        validated_input = input_schema.load(body)
    except ValidationError as e:
        print("Validation error:", e)
        return generate_response(
            400,
            {"message": e.messages}
        )

    # Now we no longer process a single SMILES. Instead we run a batch
    # from merged_db.csv and store the results.
    try:
        # filters from the request
        filters = validated_input["filters"]
        process_batch(filters)
        # After finishing, just return 200
        return generate_response(200, {"message": "Batch processing complete."})
    except Exception as e:
        print("Internal server error:", e)
        return generate_response(
            500,
            {"message": "Internal Server Error"}
        )


def generate_response(status_code, body, is_options=False):
    headers = {
        "Access-Control-Allow-Origin": "http://localhost:3001",
        "Access-Control-Allow-Methods": "GET, POST, OPTIONS",
        "Access-Control-Allow-Headers": "Content-Type, Authorization",
        "Access-Control-Allow-Credentials": "true"
    }

    if is_options:
        return {
            'statusCode': status_code,
            'headers': headers,
            'body': ''
        }

    body_str = json.dumps(body) if isinstance(body, (dict, list)) else str(body)
    return {
        'statusCode': status_code,
        'headers': headers,
        'body': body_str
    }
