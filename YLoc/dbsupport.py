import mysql.connector;


def db_query(statement):
        # create a connection to the database
        #db = MySQLdb.connect(host = "www-bs3", port = 3306, user="yloc",passwd="elucidate",db ="YLocDB")
        db = mysql.connector.connect(user='yloc', password='!yloc818%/', host='localhost', database='YLocDB')

        # get a database cursor from that connection
        cursor = db.cursor();
        cursor.execute(statement)

        if cursor.with_rows:
                return cursor.fetchall()
        else:
                return []


def getJobInfo(id):
    sql_string = "SELECT id,nr_sequences,finished_sequences FROM pending WHERE id='"+str(id)+"';";
    res=db_query(sql_string);

    if len(res) > 0:
        return res[0];
    else:
        sql_string = "SELECT id FROM queries WHERE id='"+str(id)+"';";
        res=db_query(sql_string);

        if len(res) > 0:
            return [1];
        else:
            return [];


def alterJobStatus(id, nr=1, add=True, nr_finished=-1, ip=""):
    if add:
        sql_string = "Insert Into pending (id,nr_sequences,finished_sequences,ip) values ('"+str(id)+"','"+str(nr)+"','0','"+str(ip)+"');";
        res=db_query(sql_string);
    else:
        if nr_finished > 0:
            sql_string = "Update pending set finished_sequences='"+str(nr_finished)+"' where id='"+str(id)+"';";
            res=db_query(sql_string);
        else:
            sql_string = "Delete From pending where id='"+str(id)+"';";
            res=db_query(sql_string);
