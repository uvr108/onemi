from obspy.clients.seedlink.easyseedlink import create_client

# A function to handle incoming data
def handle_data(trace):
    print('Received the following trace:')
    print(trace)
    print()

# Create the client and pass the function as a callback
client = create_client('10.54.217.11', on_data=handle_data) # SISMOLOGICAS
client.select_stream('CX', 'PB01', 'H[N/H]Z')
# client = create_client('10.54.217.11', on_data=handle_data)  # RNA
# client.select_stream('C1', 'A14C', 'HNZ')
client.run()
