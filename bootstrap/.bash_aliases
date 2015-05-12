# Define some shorthand aliases to make life easier
# @author Matthew Norris <Matthew.Norris@jic.ac.uk>

. /vagrant/bootstrap/functions.sh

# Tail apache error log
alias fa-tail='sudo tail -f /var/log/apache2/error.log'

# DL sauce
alias fa-dlSauce='( dl_sauce )'

# Hydrate database
alias fa-hydrateDB='( hydrate_db )'

# Import database
alias fa-importDB='( import_db )'

# Export database
alias fa-exportDB='( export_db )'

# Run the development server
alias fa-runDevServer='cd /vagrant/foldatlas && python3 app.py'

