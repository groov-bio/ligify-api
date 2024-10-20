FROM public.ecr.aws/lambda/python:3.12

COPY ligify/. ./

RUN python.12 -m pip install -r requirements.txt

CMD ["main.lambda_handler"]