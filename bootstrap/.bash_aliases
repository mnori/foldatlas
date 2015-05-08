# Define some shorthand aliases to make life easier
# @author Matthew Norris <Matthew.Norris@jic.ac.uk>

bash /vagrant/bootstrap/functions.sh

# Tail apache error log
alias rb-tail='sudo tail -f /var/log/apache2/error.log'

# Hydrate database
alias rb-hydrateDB='hydrate_db'

# Import database
alias rb-importDB='import_db'

# Export database
alias rb-exportDB='export_db'

# Run the development server
alias rb-runDevServer='cd /vagrant/rnabrowser && python3 app.py'

