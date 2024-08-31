from pprint import pprint
import re
from dotenv import load_dotenv
import json

from fetch_data import fetch_data
from genbank.create_genbank import create_genbank
from predict.pubchem import get_inchikey, get_name
import numpy

from marshmallow import Schema, fields, ValidationError, validate

def validate_bool(value):
    if not isinstance(value, bool):
        raise ValidationError("reviewed requires a boolean value.")

class FilterSchema(Schema):
    max_reactions = fields.Integer(required=True, validate=validate.Range(min=1, error="max_reaction must be a positive integer."))
    proteins_per_reaction = fields.Integer(required=True, validate=validate.Range(min=1, error="proteins_per_reaction must be a positive integer."))
    reviewed = fields.Boolean(required=True, validate=validate_bool)
    lineage = fields.Str(required=True, validate=validate.OneOf(['Phylum', 'Class', 'Order', 'Family', 'Genus', 'None'],
                                                                    error="lineage should be one of Domain, Phylum, Class, Order, Family, Genus, None"
                                                                ))
    max_operons = fields.Integer(required=True, validate=validate.Range(min=1, error="max_operions must be a positive integer."))
    max_alt_chems = fields.Integer(required=True, validate=validate.Range(min=1, error="max_alt_chems must be a positive integer."))

class InputSchema(Schema):
    # This might need to be adjusted in the future
    # https://gist.github.com/lsauer/1312860/63687490915641fdb4e377cf581211da9d5d64c0
    smiles = fields.String(validate=validate.Regexp(r'^([^J][a-z0-9@+\-\[\]\(\)\\\/%=#$]{6,})$', flags=re.IGNORECASE, error="Invalid SMILES provided"))
    filters = fields.Nested(FilterSchema) 

def lambda_handler(event, context):
    print(numpy.__file__)

    # def create_plasmid(regulators, chemical):
    #     for regulator in regulators:
    #         result = create_genbank(regulator['refseq'], chemical, regulator['protein']['context']['promoter']['regulated_seq'], regulator['reg_protein_seq'])
    #         regulator['plasmid_sequence'] = str(result)
        
    #     return regulators

    load_dotenv()

    body = json.loads(event['body'])

    input_schema = InputSchema()

    try:
        input_schema.load(body)
    except ValidationError as e:
        print(e)
        return {
            "status_code": 400,
            "message": json.dumps(e.messages_dict)
        }

    try:
        chemical_name = get_name(body['smiles'], 'smiles')
        InChiKey = get_inchikey(body['smiles'], 'smiles')
        chemical = {"name": chemical_name, "smiles": body['smiles'], "InChiKey": InChiKey}

        regulators, metrics = fetch_data(chemical['InChiKey'], body['filters'])
        # regulators = create_plasmid(regulators, chemical_name)
        # with open('data.json', 'w', encoding='utf-8') as f:
        #     json.dump(regulators, f, ensure_ascii=False, indent=4)
        return {
            "statusCode": 200,
            "body": json.dumps({
                # 'metrics': metrics,
                'regulators': regulators
            })
        }
    except Exception as e:
        print(e)
        # TODO - improve error responses per diagram
        return {
            "statusCode": 500,
            # "message": json.dumps()
        }

