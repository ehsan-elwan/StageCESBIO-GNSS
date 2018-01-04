#!/bin/bash

if [ $SUDO_USER ]; then uname=$SUDO_USER; else uname=`whoami`; fi
release=$(lsb_release -r)
linuxVersion=${release#*:}
versionNum=${linuxVersion%.*}

mkdir ~/maps
sudo chown $uname:$uname ~/maps
RED='\033[0;31m'
NC='\033[0m' # No Color
echo -e "Please download your map file from ${RED}http://download.geofabrik.de/${NC} "
read -n1 -r -p "to ~/maps/ folder then press space key to continue ..." key

if [ "$key" = '' ]; then

if [ $versionNum -eq 16 ]
then

apt install -y libboost-all-dev git-core tar unzip wget bzip2 build-essential autoconf libtool libxml2-dev libgeos-dev libgeos++-dev libpq-dev libbz2-dev libproj-dev munin-node munin libprotobuf-c0-dev protobuf-c-compiler libfreetype6-dev libpng12-dev libtiff5-dev libicu-dev libgdal-dev libcairo-dev libcairomm-1.0-dev apache2 apache2-dev libagg-dev liblua5.2-dev ttf-unifont lua5.1 liblua5.1-dev libgeotiff-epsg make cmake g++ libboost-dev libboost-system-dev libboost-filesystem-dev libexpat1-dev zlib1g-dev libbz2-dev libpq-dev libgeos-dev libgeos++-dev libproj-dev lua5.2 liblua5.2-dev npm nodejs-legacy php git

apt-get install -y postgresql postgresql-contrib postgis postgresql-9.5-postgis-2.2 autoconf apache2-dev libtool libxml2-dev libbz2-dev libgeos-dev libgeos++-dev libproj-dev gdal-bin libgdal1-dev libmapnik-dev mapnik-utils python-mapnik libapache2-mod-php7.0 fonts-noto-cjk fonts-noto-hinted fonts-noto-unhinted ttf-unifont python-pip


query=SQL.sql
touch $query
echo "Query Creating.."
echo "CREATE USER "$uname" ;" >> $query
echo "CREATE DATABASE gis OWNER "$uname" ENCODING 'UTF8' ;" >> $query
echo "\c gis ;" >> $query
echo "CREATE EXTENSION postgis ;" >> $query
echo "CREATE EXTENSION hstore ;" >> $query
echo "ALTER TABLE geometry_columns OWNER TO "$uname" ;" >> $query
echo "ALTER TABLE spatial_ref_sys OWNER TO "$uname" ;" >> $query
sudo -u postgres bash -c "psql -f SQL.sql"
rm $query


mkdir ~/src
sudo chown $uname:$uname ~/src
cd ~/src
git clone git://github.com/openstreetmap/osm2pgsql.git
cd osm2pgsql
mkdir build && cd build
cmake ..
make
sudo make install

cd ~/src
git clone git://github.com/SomeoneElseOSM/mod_tile.git
cd mod_tile
./autogen.sh
./configure
make
sudo make install
sudo make install-mod_tile
sudo ldconfig

cd ~/src
git clone git://github.com/gravitystorm/openstreetmap-carto.git
cd openstreetmap-carto
npm install -g carto
carto project.mml > mapnik.xml

key=1

cd ~/maps

file=$(ls *.pbf)

su $uname -c "osm2pgsql -d gis --create --slim  -G --hstore --tag-transform-script ~/src/openstreetmap-carto/openstreetmap-carto.lua -C 2500 --number-processes 4 -S ~/src/openstreetmap-carto/openstreetmap-carto.style $file"


cd ~/src/openstreetmap-carto/

scripts/get-shapefiles.py

path="XML=/home/"$uname"/src/openstreetmap-carto/mapnik.xml"
old="XML=/home/renderaccount/src/openstreetmap-carto/mapnik.xml"
sudo sed -i 's,'$old','$path',g' /usr/local/etc/renderd.conf


sudo mkdir /var/lib/mod_tile
sudo chown $uname /var/lib/mod_tile
sudo mkdir /var/run/renderd
sudo chown $uname /var/run/renderd


apachconf="/etc/apache2/sites-available/000-default.conf"
sudo touch /etc/apache2/conf-available/mod_tile.conf
echo "LoadModule tile_module /usr/lib/apache2/modules/mod_tile.so" >> /etc/apache2/conf-available/mod_tile.conf
sudo a2enconf mod_tile
sudo sed -i '12i    ModTileMissingRequestTimeout 30' $apachconf
sudo sed -i '12i    # Timeout before giving up for a tile to be rendered that is otherwise missing' $apachconf
sudo sed -i '12i    ModTileRequestTimeout 0' $apachconf
sudo sed -i '12i    # Timeout before giving up for a tile to be rendered' $apachconf
sudo sed -i '12i    ModTileRenderdSocketName /var/run/renderd/renderd.sock' $apachconf
sudo sed -i '12i    LoadTileConfigFile /usr/local/etc/renderd.conf' $apachconf


sudo sed -i '19i</Directory>' $apachconf
sudo sed -i '19iOptions +ExecCGI' $apachconf
sudo sed -i '19i<Directory /home/'$uname'/Spec_interface/python/>' $apachconf
sudo sed -i '19i</Directory>' $apachconf
sudo sed -i '19iAddHandler cgi-script .py' $apachconf
sudo sed -i '19iallow from all' $apachconf
sudo sed -i '19iOrder allow,deny' $apachconf
sudo sed -i '19iAllowOverride None' $apachconf
sudo sed -i '19iOptions ExecCGI Indexes FollowSymLinks MultiViews' $apachconf
sudo sed -i '19i<Directory /var/www/>' $apachconf

sudo service apache2 reload

path="RUNASUSER="$uname
sudo sed -i 's,RUNASUSER=renderaccount,'$path',g' ~/src/mod_tile/debian/renderd.init
sudo cp ~/src/mod_tile/debian/renderd.init /etc/init.d/renderd
sudo chmod u+x /etc/init.d/renderd
sudo cp ~/src/mod_tile/debian/renderd.service /lib/systemd/system/
sudo /etc/init.d/renderd start
sudo systemctl enable renderd

sudo pip install numpy
sudo pip install pandas
sudo pip install geojson
sudo pip install requests
sudo pip install pynmea2
sudo pip install openpyxl
sudo pip install ephem
sudo pip install pyserial
sudo pip install shapely
sudo pip install geopandas
sudo pip install pyproj

cd ~/src
git clone https://github.com/tkrajina/srtm.py.git
cd srtm.py
python setup.py build
sudo python setup.py install

cd ~
sudo mv stage-ehsan Spec_interface
sudo chown -R $uname:$uname Spec_interface

path='"/home/'$uname'/Spec_interface/python/SRTM_cache/"'
main="/usr/local/lib/python2.7/dist-packages/srtm/main.py"
var='return '$path

sudo sed -i 's,result = "",result = '"$path"',g' $main
sudo sed -i 's,return result,return '"$path"',g' $main
sudo sed -i "s,raise Exception('No,pass #raise Exception('No,g" $main

sudo cp -a ~/Spec_interface/html/. /var/www/html/
sudo ln -s ~/Spec_interface/ /var/www/html/
sudo ln -s ~/Spec_interface/python/ /var/www/html/
sudo chown -R www-data:www-data /var/www/html/
sudo chmod -R 777 ~/Spec_interface
sudo a2enmod cgid
cd /etc/apache2/mods-enabled
sudo ln -s /etc/apache2/mods-available/cgi.load
sudo service apache2 restart

sudo rm /var/crash/*

echo -e " ${RED}Done! enjoy....${NC} "

else
  echo "17 Version!"
fi

fi
