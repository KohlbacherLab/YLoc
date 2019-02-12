python /webservice/ylsetup.py

chown -R mysql:mysql /var/lib/mysql
/etc/init.d/mysql start

mysql -e "CREATE DATABASE YLocDB"
mysql YLocDB < /webservice/ylocdb.sql
mysql -e "CREATE USER yloc@localhost IDENTIFIED BY '!yloc818%/';"
mysql -e "GRANT ALL PRIVILEGES ON YLocDB.* TO 'yloc'@'localhost';"

/etc/init.d/apache2 start

tail -f /var/log/apache2/error.log
