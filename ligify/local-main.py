from pprint import pprint
from dotenv import load_dotenv
import json
import csv
import os
import boto3
from boto3.dynamodb.conditions import Key

from fetch_data import fetch_data
from genbank.create_genbank import create_genbank
from predict.pubchem import get_inchikey, get_name

def ensure_table_exists(dynamodb):
    try:
        table = dynamodb.Table("Chemicals")
        table.table_status
        print("Table exists:", table.table_name)
        return table
    except Exception as e:
        print("Table doesn't exist, creating...")
        try:
            table = dynamodb.create_table(
                TableName='Chemicals',
                KeySchema=[
                    {
                        'AttributeName': 'SMILES',
                        'KeyType': 'HASH'
                    },
                    {
                        'AttributeName': 'chunk_index',
                        'KeyType': 'RANGE'
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
            table.meta.client.get_waiter('table_exists').wait(TableName='Chemicals')
            print("Table created successfully!")
            return table
        except Exception as create_error:
            print("Error creating table:", create_error)
            raise

def store_in_dynamodb(table, smiles, data, data_is_error=False):
# Check if entry already exists
    try:
        existing = table.query(
            KeyConditionExpression=Key('SMILES').eq(smiles)
        )
        if existing['Items']:
            print(f"Entry already exists for SMILES: {smiles}")
            return
    except Exception as e:
        print(f"Error checking for existing entry: {e}")    

    data_str = json.dumps(data)
    data_size = len(data_str.encode('utf-8'))

    max_chunk_size = 350000
    if data_size <= max_chunk_size:
        table.put_item(
            Item={
                "SMILES": smiles,
                "chunk_index": 0,
                "data": data if not data_is_error else None,
                "error": data if data_is_error else None
            }
        )
    else:
        chunks = []
        start = 0
        while start < data_size:
            end = start + max_chunk_size
            chunk_data = data_str[start:end]
            chunks.append(chunk_data)
            start = end
        
        for i, chunk in enumerate(chunks):
            table.put_item(
                Item={
                    "SMILES": smiles,
                    "chunk_index": i,
                    "data": chunk
                }
            )

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

def process_batch(table, filters):
    csv_path = "./merge_db.csv"
    if not os.path.exists(csv_path):
        print("merge_db.csv not found!")
        return

    results = []
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        entries = list(reader)[:1000]  # first 100 entries
        for index, entry in enumerate(entries):
            print(f"current index: {index + 1}")
            smiles = entry["SMILES"]
            chemical_name = entry["Name"]
            chebid = entry["CHEBID"]
            
            print(f"Processing: {smiles}")

            try:
                try:
                    inchi_key = get_inchikey(smiles, "smiles")
                except Exception as e:
                    error_item = {
                        "error": f"Failed to get InChIKey: {str(e)}"
                    }
                    store_in_dynamodb(table, smiles, error_item, True)
                    continue

                try:
                    regulators, metrics = fetch_data(inchi_key, filters)
                except Exception as e:
                    error_item = {
                        "error": f"Failed to fetch data: {str(e)}"
                    }
                    store_in_dynamodb(table, smiles, error_item, True)
                    continue

                try:
                    regulators = create_plasmid(regulators, chemical_name)
                except Exception as e:
                    error_item = {
                        "error": f"Failed to create plasmid: {str(e)}"
                    }
                    store_in_dynamodb(table, smiles, error_item, True)
                    continue

                full_item = {
                    "CHEBID": chebid,
                    "SMILES": smiles,
                    "Name": chemical_name,
                    "metrics": metrics,
                    "regulators": regulators
                }
                store_in_dynamodb(table, smiles, full_item)

            except Exception as e:
                error_item = {
                    "SMILES": smiles,
                    "error": f"Unexpected error: {str(e)}"
                }
                store_in_dynamodb(table, smiles, error_item, True)

    return results

if __name__ == "__main__":
    load_dotenv()
    
    # DynamoDB setup
    dynamodb = boto3.resource(
        'dynamodb',
        endpoint_url="http://localhost:8000",
        region_name="localhost",
        aws_access_key_id="fakeMyKeyId",
        aws_secret_access_key="fakeSecretAccessKey"
    )
    
    table = ensure_table_exists(dynamodb)
    
    # Define filters
    filters = {
        "max_reactions": 20,
        "proteins_per_reaction": 20,
        "reviewed": True,
        "lineage": 'None',
        "max_operons": 20,
        "max_alt_chems": 10
    }

    try:
        results = process_batch(table, filters)
        
        # Save results to a JSON file
        # with open('batch_results.json', 'w', encoding='utf-8') as f:
        #     json.dump(results, f, ensure_ascii=False, indent=4)
        
        print("Batch processing completed. Results saved to batch_results.json")
        
    except Exception as e:
        print(f"Error during batch processing: {e}")