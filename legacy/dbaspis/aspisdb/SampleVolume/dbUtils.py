import mysql.connector
from mysql.connector import Error
import configparser

class dbAccess():
 
 def __init__(self):
  self.config = {}
  self.hasTransaction = False
  self.lastError = ""
  self.connection = None
  
#-------------------------------------------------
 def loadConfig(self,path = "con.fig"):
  c = configparser.ConfigParser()
  c.read(path)
  for section in c.sections():
   print(section)
   if section == "mysql":
    for key in c[section]:
     if key != "password":
      print("{} = [{}]".format(key,c[section][key]))
     else:
      print("{} = [{}]".format(key,"x" * 6))
     self.config[key] = c[section][key]
#-------------------------------------------------
 def getDbConnection(self,_keepConnection = False):
  
  if self.hasTransaction:
   return self.connection
  
  try:
   
   if _keepConnection:
    self.connection = mysql.connector.connect(user=self.config['user'],
                                  host=self.config['host'],
                                  port=self.config['port'],
                                  password=self.config['password'],
                                  database=self.config['database'])
    return self.connection
   else:
    return mysql.connector.connect(user=self.config['user'],
                                  host=self.config['host'],
                                  port=self.config['port'],
                                  password=self.config['password'],
                                  database=self.config['database'])

  except Error as e:
   print("Failed to connect to database")
   print(e)
   return None

#-------------------------------------------------  
 def query(self, sql, args):
  conn = self.getDbConnection()
  if conn:
   try:
    cursor = conn.cursor()
    cursor.execute(sql,args)
    results = cursor.fetchall()
    cursor.close()
    if not self.hasTransaction:
     conn.close()
    return results
   except Exception as ex:
    self.lastError = ex
    print("DB ERROR: {}".format(ex))
    return None
  else:
   return None
   
#-------------------------------------------------  
 def execute(self, sql, args):
  conn = self.getDbConnection()
  if conn:
   try:
    cursor = conn.cursor()
    cursor.execute(sql,args,multi=True)
    insertId = cursor.lastrowid
    
    cursor.close()
    if not self.hasTransaction:
     conn.commit()
     conn.close()
     
    return (True,insertId)
   except Exception as ex:
    self.lastError = ex
    return (False,"Error")
  else:
   return (False,"No connection")
   
#------------------------------------------------- 
 def beginTransaction(self):
  if self.hasTransaction:
   return (False,"A transaction has already been started")
   
  try:
   self.getDbConnection(True)
   self.hasTransaction = True
   return (True, "transaction started")
  except Exception:
   return (False,"Failed to start transaction")
   
#------------------------------------------------- 
 def commitTransaction(self):
  if not self.hasTransaction:
   return (False,"A transaction has not been started")
   
  try:
   conn = self.getDbConnection()
   conn.commit()
   self.hasTransaction = False
   conn.close()
   return (True, "Commited")
  except Exception:
   return (False,"Failed to commit")
   
#------------------------------------------------- 
 def rollbackTransaction(self):
  if not self.hasTransaction:
   return (False,"A transaction has not been started")
   
  try:
   conn = self.getDbConnection()
   conn.rollback()
   self.hasTransaction = False
   conn.close()
   return (True, "rolled back")
  except Exception:
   return (False,"Failed to rollback")   
 
#------------------------------------------------- 
 def __del__(self):
  if self.hasTransaction:
   print("ERROR! Transaction has not been commited! Rolling back")
   self.rollbackTransaction()
  
#------------------------------------------------- 
 def makeSQLCompatibleName(self,_name):
  return _name.replace(":","_").replace(" ","_").replace(".","_").replace("-","_")
   
