import re
from dotenv import load_dotenv
import json

from fetch_data import fetch_data
from genbank.create_genbank import create_genbank
from predict.pubchem import get_inchikey, get_name

from marshmallow import Schema, fields, ValidationError, validate


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


def lambda_handler(event, context):
    print("Received event:", json.dumps(event))

    # Extract path and method from the event
    path = event.get('rawPath') or event.get('path')
    method = event.get('httpMethod', '').upper()
    print(f"Method: {method}, Path: {path}")  # Log the method and path

    # Adjust path validation
    if not path or '/ligify' not in path:
        return generate_response(403, "Forbidden")

    # Handle CORS preflight request
    if method == 'OPTIONS':
        return generate_response(200, "", is_options=True)

    # Load environment variables
    load_dotenv()

    # Parse and validate the request body
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

    # Process the request
    try:
        chemical_name = get_name(validated_input["smiles"], "smiles")
        InChiKey = get_inchikey(validated_input["smiles"], "smiles")
        chemical = {
            "name": chemical_name,
            "smiles": validated_input["smiles"],
            "InChiKey": InChiKey,
        }

        regulators, metrics = fetch_data(chemical["InChiKey"], validated_input["filters"])
        regulators = create_plasmid(regulators, chemical_name)

        response_body = {
            "metrics": metrics,
            "regulators": regulators
        }

        return generate_response(200, response_body)
    except Exception as e:
        print("Internal server error:", e)
        return generate_response(
            500,
            {"message": "Internal Server Error"}
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


def generate_response(status_code, body, is_options=False):
    headers = {
        "Access-Control-Allow-Origin": "http://localhost:3001",
        "Access-Control-Allow-Methods": "GET, POST, OPTIONS",
        "Access-Control-Allow-Headers": "Content-Type, Authorization",
        "Access-Control-Allow-Credentials": "true"
    }

    if is_options:
        # OPTIONS requests typically do not have a body
        return {
            'statusCode': status_code,
            'headers': headers,
            'body': ''
        }

    # Ensure the body is a JSON string
    body_str = json.dumps(body) if isinstance(body, (dict, list)) else str(body)

    return {
        'statusCode': status_code,
        'headers': headers,
        'body': body_str
    }
