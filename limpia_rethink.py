from rethinkdb import RethinkDB

r = RethinkDB()
rethink_r = r.connect(host='10.54.223.19', port=28015, db='csn').repl()


res = r.table('ULI_SOFT').delete().run()

print(res)

rethink_r.close()
