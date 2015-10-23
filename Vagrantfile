# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure(2) do |config|
    config.vm.box = "ubuntu/trusty64"
    config.vm.provision :shell, path: "bootstrap/bootstrap.sh"
    config.vm.network "private_network", ip: "192.168.50.2"
    config.vm.provider "virtualbox" do |v|
        v.customize ["modifyvm", :id, "--memory", 4096, "--ioapic", "on", "--cpus", 2]

        # this sets up a shared folder in /vagrant/sauce_data, pointing to all the source data. useful for importing stuff
        config.vm.synced_folder "~/data/foldatlas", "/vagrant/sauce_data"

        # JIC setup for structures
        # another folder will point to the structures
        
        # config.vm.synced_folder "/media/shares/Research-Groups/Yiliang-Ding/data_analysis_Ding_2013/MAC/Yin/Mapping_F/raw_data/structures", "/vagrant/structure_data"

        # home setup for structures
        config.vm.synced_folder "~/data_input/foldatlas_structures", "/vagrant/structure_data"
    end
end
