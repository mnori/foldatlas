# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure(2) do |config|
    config.vm.box = "ubuntu/trusty64"
    config.vm.provision :shell, path: "bootstrap/bootstrap.sh"
    config.vm.network "private_network", ip: "192.168.50.2"
    config.vm.provider "virtualbox" do |v|
        v.customize ["modifyvm", :id, "--memory", 4096, "--ioapic", "on", "--cpus", 2]
    end
end
