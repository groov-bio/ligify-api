FROM public.ecr.aws/lambda/python:3.11

COPY ligify/. ./

RUN python3.11 -m pip install -r requirements.txt

CMD ["main.lambda_handler"]