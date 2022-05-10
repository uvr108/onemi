from obspy.clients.seedlink.easyseedlink import create_client

# A function to handle incoming data
def handle_data(trace):
    print('Received the following trace:')
    print(trace)
    print()

# Create the client and pass the function as a callback
client = create_client('10.54.217.12', on_data=handle_data)
client.select_stream('CX', 'PB01', 'HHZ')
client.run()
