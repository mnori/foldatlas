# VERY USEFUL
# https://gist.github.com/jmorton/558f7079ed2159156277

# This kind of thing will work:

# This line is key - points to the Python3.4 WSGI apache module
# If you load the wrong one, you have to remove the conf and the .so module itself.
LoadModule wsgi_module /etc/httpd/modules/mod_wsgi-py34.cpython-34m.so

# Main foldatlas
<VirtualHost *>
    ServerName v0622.nbi.ac.uk
    WSGIDaemonProcess foldatlas user=norrism group=norrism threads=5
    
    # \ python-path=/usr/local/lib/python3.5/dist-packages
                
    # Serve foldatlas from a subfolder off the main domain
    WSGIScriptAlias / /var/www/foldatlas/foldatlas/foldatlas.wsgi

    <Directory /var/www/foldatlas/foldatlas/>
        WSGIProcessGroup foldatlas
        WSGIApplicationGroup %{GLOBAL}
        Require all granted
    </Directory>
</VirtualHost>

# no need to map static - it's in the parent folder already :D

# vim: sntax=apache ts=4 sw=4 sts=4 sr noet

foldatlas_wsgi.conf (END)



/usr/local/bin/pip3.4

HOST: v0622.nbi.ac.uk
USER: norrism
PASS: 5Y5THNR2

sudo su

# install apache
sudo yum install httpd
sudo service httpd start

# now, add these lines to apache config
LoadModule proxy_module modules/mod_proxy.so
LoadModule proxy_http_module modules/mod_proxy_http.so
LoadModule headers_module modules/headers.so

# install mysql
sudo yum install mysql-server
sudo service mysqld start

# set the password
mysql -uroot
SET PASSWORD FOR 'root'@'localhost' = PASSWORD('jGEHL3qT6sdntJD9pfyB8f3hGzBajLW2');
create database foldatlas;
echo "create database foldatlas" | mysql -u root -p

# sudo yum install phpmyadmin # nope.jpeg

# Install python3
https://www.softwarecollections.org/en/scls/rhscl/rh-python34/

# run this before literally everything - how to make automatic
scl enable rh-python34 bash

# sklearn
pip3 install numpy
pip3 install scipy
conda install scikit-learn




