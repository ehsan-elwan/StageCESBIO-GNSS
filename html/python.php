
<?php
$target_dir = "Spec_interface/GeoJson/";
$target_file = $target_dir . "Obslayer.json";
$uploadOk = 1;
$FileType = pathinfo($target_file,PATHINFO_EXTENSION);
// Check if image file is a actual image or fake image

// Allow certain file formats
if($FileType != "geojson" && $FileType != "json" ) {
    echo "Sorry, only GeoJSON and JSON files are allowed.";
    $uploadOk = 0;
}

if (file_exists($target_file)) {
    ; //Change the file permissions if allowed
    unlink($target_file); //remove the file
}
// Check if $uploadOk is set to 0 by an error
if ($uploadOk == 0) {
    echo "Sorry, your file was not uploaded.";
// if everything is ok, try to upload file
} else {
    if (move_uploaded_file($_FILES["FileInput"]["tmp_name"], $target_file)) {
        echo "The file ". basename( $_FILES["FileInput"]["name"]). " has been uploaded.";
        chmod($target_file,0777);
        exec('python Spec_interface/python/Obs_reader.py');
        echo "passed";
    } else {
        echo "Sorry, there was an error uploading your file.";
    }
}
?>
