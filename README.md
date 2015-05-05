# RNA Browser
Matthew Norris <Matthew.Norris@jic.ac.uk>

## Prerequisites

- Vagrant 1.7.2 or above
- Preferably OSX or Linux 

## Setting up the box

Run `vagrant up` in the same folder as this file. The first `vagrant up` call will install everything and download genome sequence data. If it fails, first check your `vagrant --version`. Failing that, pester the author.

## Hosts

	192.168.50.2	rnabrowser.dev
	192.168.50.2    static.rnabrowser.dev
	192.168.50.2    pma.rnabrowser.dev

Add these lines to your `hosts` file, typically located in /etc/hosts on *nix based systems.

## Running the web server

    vagrant ssh
    rb-runDevServer

This runs a Flask development web server. Go to http://rnabrowser.dev to see the site running. The server will keep going until you press `CTRL-C` in the terminal.

## Resetting the database

    vagrant ssh
    rb-resetDB

This drops all the database tables, recreates them, and then populates the tables by parsing the *.fa and *.gff3 files downloaded during the `vagrant up` provisioning. There will soon be lots of delicious metadata attached to each transcript, describing RNA structure predictions, experimental data, etc.

## Viewing / hacking the database

Go to http://pma.rnabrowser.dev

    username: root
    password: vagrant


