# Define some shorthand aliases to make life easier
# @author Matthew Norris <Matthew.Norris@jic.ac.uk>

# Tail apache error log
alias rb-tail='sudo tail -f /var/log/apache2/error.log'

# Reset database
alias rb-resetDB='python3 /vagrant/rnabrowser/app.py resetdb '

# Run the development server
alias rb-runDevServer='cd /vagrant/rnabrowser && python3 app.py'
