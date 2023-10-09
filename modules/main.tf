
//TODO zone region subnetwork　machine_typeの外だし
provider "google" {
  credentials = file("life-science-1st-project.json")

  project = "life-science-1st-project"
  region  = var.region
  zone    = var.zone
}

resource "google_compute_instance" "dev-avocado" {
  name = format("%s-%s-maintantest-terraform-instance%02d", var.project_name, var.case, count.index)
  #custom machine type 6 for vcpu core and 20480 /1024 = 20GB memory
  #limit 6.5 GB per CPU unless you add extended memory. add suffix -ext for extra memory i.e. custom-2-15360-ext
  #https://cloud.google.com/dataproc/docs/concepts/compute/supported-machine-types
  #machine_type = "${format("custom-%s-%s", local.instance_cpu, local.instance_memory)}"
  machine_type = var.machine_type
  #count        = tonumber("${var.seeds}") / tonumber(local.instance_cpu) / 10
  count = tonumber("${var.seeds}") / var.machine_ratio
  tags  = ["externalssh"]
  boot_disk {
    mode        = "READ_WRITE"
    auto_delete = true
    device_name = "dev-avocado-test"
    initialize_params {
      size = 100
      #  image = "image-combo-1"
      #2022-10-17 dev-avacoda logo
      image = "physbo-tempalte-0227"
    }
  }


  network_interface {
    subnetwork = var.subnetwork
    #  access_config {
    #    nat_ip = google_compute_address.static[count.index].address
    #  }
  }
  scheduling {
    on_host_maintenance = "MIGRATE"
    provisioning_model  = "STANDARD"
  }

  service_account {
    scopes = ["https://www.googleapis.com/auth/servicecontrol", "https://www.googleapis.com/auth/service.management.readonly", "https://www.googleapis.com/auth/logging.write", "https://www.googleapis.com/auth/monitoring.write", "https://www.googleapis.com/auth/trace.append", "https://www.googleapis.com/auth/devstorage.full_control"]
  }

  shielded_instance_config {
    enable_secure_boot          = false
    enable_vtpm                 = true
    enable_integrity_monitoring = true

  }

  reservation_affinity {
    type = "ANY_RESERVATION"
  }


  #metadata_startup_script = "${format("/usr/bin/python3.7 /home/atom/work/revolka/combo/rk002_ml_test/script/remote_run.py %02d %s", count.index, var.project_name)}"
  #metadata_startup_script =file("runtest.sh")
  #make sure no extra space

  metadata_startup_script = file("run.sh")
  metadata = {
    enable-oslogin = true
   # ssh-keys       = "atom:${file("id_rsa.pub")}"
  }


}

output "debug_print" {
  value = "${var.case} to ${var.project_name}"
}
#~ ❯❯❯ for i in `echo {a..b}`
#do
#gcloud compute instances create dev-avocado-$i --project=life-science-1st-project 
#--zone=asia-northeast1-b --machine-type=c2-standard-60 --network-interface=subnet=tokyo,no-address --metadata=enable-oslogin=true --maintenance-policy=MIGRATE --provisioning-model=STANDARD --service-account=927017471470-compute@developer.gserviceaccount.com --scopes=https://www.googleapis.com/auth/servicecontrol,https://www.googleapis.com/auth/service.management.readonly,https://www.googleapis.com/auth/logging.write,https://www.googleapis.com/auth/monitoring.write,https://www.googleapis.com/auth/trace.append,https://www.googleapis.com/auth/devstorage.full_control --create-disk=auto-delete=yes,boot=yes,device-name=dev-avocado-$i,image=projects/life-science-1st-project/global/images/dev-avocado-im-220515,mode=rw,size=100,type=projects/life-science-1st-project/zones/asia-northeast1-b/diskTypes/pd-balanced 
#--no-shielded-secure-boot --shielded-vtpm --shielded-integrity-monitoring --reservation-affinity=any
#done
