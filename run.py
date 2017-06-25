#!flask/bin/python
from app import app
from app.edit import init_data_vars
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

init_data_vars()
if __name__ == '__main__':
    #deployment definitions. 
    app.run(host='0.0.0.0', port=40001, debug=True)
    #for local deployment. for tests
    #app.run(host='127.0.0.1', port=5000, debug=True)
