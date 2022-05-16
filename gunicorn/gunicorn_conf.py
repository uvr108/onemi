from multiprocessing import cpu_count

# Socket Path
bind = 'unix:/home/ulises/Documentos/Projectos/Python/onemi/gunicorn/gunicorn.sock'

# Worker Options
workers = cpu_count() + 1
worker_class = 'uvicorn.workers.UvicornWorker'

# Logging Options
loglevel = 'debug'
accesslog = '/home/ulises/Documentos/Projectos/Python/onemi/gunicorn/access_log'
errorlog =  '/home/ulises/Documentos/Projectos/Python/onemi/gunicorn/error_log'
