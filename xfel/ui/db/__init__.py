from __future__ import division
try:
  import MySQLdb
except ImportError, e:
  from libtbx.utils import Sorry
  raise Sorry("Mysql is not installed")

def get_run_path(rootpath, trial, rungroup, run):
  import os
  return os.path.join(rootpath, "r%04d"%run.run, "%03d_rg%03d"%(trial.trial, rungroup.id))

def get_db_connection(params, block=True):
  if params.db.password is None:
    password = ""
  else:
    password = params.db.password
  while True:
    try:
      dbobj=MySQLdb.connect(passwd=password,user=params.db.user,host=params.db.host,db=params.db.name)
      return dbobj
    except Exception,e:
      if "Too many connections" in e and block:
        print "Too many connections.  Blocking..."
        dbobj = None
        import time
        time.sleep(1)
        continue
      else:
        raise e

class db_proxy(object):
  def __init__(self, app, table_name, id = None, **kwargs):
    self.db_dict = {}
    self.app = app
    self.id = id
    self.table_name = table_name

    if id is None:
      # add the new items to the db
      query = "INSERT INTO `%s` " % self.table_name
      keys = []
      vals = []
      for key, value in kwargs.iteritems():
        self.db_dict[key] = value
        keys.append(key)
        if isinstance(value, bool):
          value = "'%s'"%int(value)
        elif value is None or value == "None" or value == "":
          value = "NULL"
        else:
          value = "'%s'"%str(value)
        vals.append(value)
      query += "(%s) VALUES (%s)"%(", ".join(keys), ", ".join(vals))
      cursor = self.app.execute_query(query, commit=True)
      self.id = cursor.lastrowid
    else:
      query = "SHOW COLUMNS FROM `%s`" % self.table_name
      cursor = self.app.execute_query(query)
      columns = [c[0] for c in cursor.fetchall()]

      query = "SELECT * FROM `%s` WHERE id = %d" % (self.table_name, id)
      cursor = self.app.execute_query(query)
      data = cursor.fetchall()[0]

      for key, value in zip(columns, data):
        if key == 'id':
          continue
        self.db_dict[key] = value

  def __getattr__(self, key):
    # Called if the property is not found
    if key in self.db_dict:
      return self.db_dict[key]
    else:
      raise AttributeError()

  def __setattr__(self, key, value):
    # Need to test for db_dict to avoid infinite loop
    #  Test key == "db_dict" to allow creating db_dict
    #  Test hasattr(self, "db_dict") in case dbproxy.__init__ hasn't been called
    #  Test key not in self.db_dict to allow setting member variables not in database dictionary
    if key == "db_dict" or not hasattr(self, "db_dict") or key not in self.db_dict:
      super(db_proxy, self).__setattr__(key, value)
      return

    self.db_dict[key] = value

    if isinstance(value, bool):
      value = "%s"%int(value)
    elif value is None or value == "None" or value == "":
      value = "NULL"
    else:
      value = "'%s'"%value

    query = "UPDATE `%s` SET %s = %s WHERE id = %d"% (
      self.table_name, key, value, self.id)
    print query
    self.app.execute_query(query, commit=True)
