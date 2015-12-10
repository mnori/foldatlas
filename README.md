# FoldAtlas
Matthew Norris <Matthew.Norris@jic.ac.uk>

## Prerequisites

- Latest version of Vagrant and Virtualbox
- Preferably OSX or Linux 

## Setting up the box

Run `vagrant up` in the same folder as this file. The first `vagrant up` call will install everything and download genome sequence data. If it fails, first check your `vagrant --version`.

## Hosts

	192.168.50.2	foldatlas.dev
	192.168.50.2    static.foldatlas.dev
	192.168.50.2    pma.foldatlas.dev

Add these lines to your `hosts` file, typically located in /etc/hosts on *nix based systems.

## Running the web server

    vagrant ssh
    rb-runDevServer

This runs a Flask development web server. Go to http://foldatlas.dev to see the site running. The server will keep going until you press `CTRL-C` in the terminal.

## Resetting the database

    vagrant ssh
    rb-resetDB

This drops all the database tables, recreates them, and then populates the tables by parsing the *.fa and *.gff3 files downloaded during the `vagrant up` provisioning. 

## Viewing / hacking the database

Go to http://pma.foldatlas.dev

    username: root
    password: vagrant


