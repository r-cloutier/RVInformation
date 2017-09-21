<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<link rel="stylesheet" type="text/css" href="style.css" />
<title>SOSS Simu1D </title>

<?php
	$basefolder = basename(getcwd())."/";
	$baseurl = "http://maestria.astro.umontreal.ca/niriss/";
	$url = $baseurl.$basefolder;
?>

<!--
<script type="text/javascript">

  var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-34408141-1']);
  _gaq.push(['_trackPageview']);

  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();

</script>
-->
</head>

<body>

<div id="container">
	<div id="mainpic">
        	<h1>SOSS Simulation Tools</h1>
			<h2>by the NIRISS Instrument Team</h2>
        </div>   
        
        <div id="menu">
        	<ul>
            	<li class="menuitem"><a href="<?php echo $url?>simu1D.php">SOSS 1D<br>Simulator</a></li>
                <li class="menuitem"><a href="<?php echo $url?>simu2D.php">SOSS 2D<br>Simulator</a></li>
                <li class="menuitem"><a href="<?php echo $url?>../SOSS_cont/SOSScontam.php">SOSS Trace<br>Contamination</a></li>
				<li class="menuitem"><a href="<?php echo $url?>../AMIcontrast/index.php">AMI Contrast<br><a></li>
            </ul>
        </div>
        

		<div id="content">
			<h2>The NIRISS Single Object Slitless Spectrograph (SOSS) 1-D Simulator: </h2> 
			<h1>This tool prepared by the NIRISS Instrument Team simulates spectra produced with the SOSS observing mode of NIRISS doing transit spectroscopy. It computes signal and noise
			from first principles, based on our best knowledge of the instrument and observatory. This is a 1-D tool, i.e., we deal with flux per pixel row extracted from the detector and
			assume that the extraction was done optimally, not introducing errors. We encourage the user to read the Simulator Guide. </h1> 
			<h3>(Version 1.0)</h3> 
		</div>
		
		<div id="p">
		<h2>Here is a <a href="<?php echo $url?>SOSS_Simulator_Guide.pdf">SOSS Simulator Guide</a> (in pdf).</h2>

<!--
<h2><font color="#FF0000">WARNING! MODIFICATIONS TO THE SIMULATOR ARE BEING TESTED AT THE MOMENT. BE BACK IN A FEW MINUTES.</font></h2>
-->



					
	<?php
		//piece of code that populates an array of plnet atmosphere models
		$directory = getcwd()."/planet_models/";
		$modellist = array();//everything in the directory that is not . or ..
		foreach(glob($directory."*.fits") as $model) {
			$modellist[] = basename( $model);
		}
		$nmodel = count($modellist);

		//piece of code that presets a few simulation input parameters to defualt values
		$out_to_in_factor = 1;
		$pixelbin_m1 = 2;	//nyquist
		$pixelbin_m2 = 2;	//nyquist
		$ntransit = 1;
		//$scale_atmosphere = 1;
		$nread = 2;			//let the script find optimal nread CHANGE AGAIN TO 0 ONCE THE CLEVER ALGO EXISTS
		$readoutmode = 256;	//standard readout mode
	?>
	<?php
		if(isset($_GET['Resolve'])){
			$Data_get = shell_exec(getcwd().'/resolverphp.py'.' "'.$_GET['targetname'].'" ');
			$Results = json_decode($Data_get, true);
			//I'M SORRY THIS IS UGLY
			if($Results !== NULL ){
				$star_radius = $Results[4]["st_rad"];
				$star_teff = $Results[3]["st_teff"];
				$star_jvega = $Results[7]["st_j"];
				$transit_duration = $Results[6]["pl_trandur"];
				$planet_density = $Results[2]["pl_dens"];
				$planet_teff = $Results[5]["pl_eqt"];
				$solid_planet_radius = $Results[1]["pl_radj"];
				$spr_unit = "Rjup";
			} else {
				$star_radius = "" ;
				$star_teff = "";
				$star_jvega = "";
				$transit_duration = "";
				$planet_density = "";
				$planet_teff =	"";
				$solid_planet_radius = "";
				$spr_unit = "Rjup";
			}
		}
	?>
<div id="p">

<form action="<?php echo $url?>simu1D.php" method="get" autocomplete="off">
	<table border="1">
		<tr><td colspan=2><b>Astrophysical Inputs</b></td></tr>
		<tr><td>Planet Name</td>
		<td><input type="text" name="targetname" value="<?php echo isset($_GET['targetname']) ? $_GET['targetname'] : $targetname ?>"  size="18" maxlength="30"/>
		<input type = submit value = "Press to Resolve" name = "Resolve">
		<?php 
			if(isset($_GET['Resolve'])) { 
				if($Results == NULL) { 
					echo "failed to resolve";
				}
			}
		?>
		<div id = "p">
		<FONT COLOR ="990000"><b><u>NOTE</u>: To resolve, the same formatting as the <a href = "http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&select=pl_name&where=pl_tranflag=1">NASA Exoplanet
		Archive</a> might be needed (ex : CoRoT-1 b).</b></FONT></div></td></tr>
		<?php
			$formPlanetModel_selected = $_GET['formPlanetModel']; 
		?>
		<tr><td>Planet Atmosphere Model</td>
		<td><select name="formPlanetModel">
			<option value="HD189733b_fortney2014.fits">Generic Hot Jupiter (HD189733b_fortney2014.fits)</option>
			<option value="Earth_kaltenegger2009.fits">Generic Earth (Earth_kaltenegger2009.fits)</option>
			<option value="GJ1214b_H2O_millerricci2010.fits">Generic GJ1214b (GJ1214b_H2O_millerricci2010.fits)</option>
			<option value="flatspectrum.fits">Flat Spectrum (flatspectrum.fits)</option>
			<option value="">--------------------</option>
			<?php
				$n = 0;
				while($n < ($nmodel)){
					//echo $formPlanetModel_selected;
					if ($formPlanetModel_selected == $modellist[$n]){
			?>
						<option selected value="<?php echo $modellist[$n]; ?>"><?php echo $modellist[$n]; ?></option> 
			<?php 
					} else { 
			?> 
						<option value="<?php echo $modellist[$n]; ?>"><?php echo $modellist[$n]; ?></option>
			<?php
					}
					$n++;
				}
			?>

		</select></td>
		</tr>
		<tr>
		<td>&nbsp; </td><td><input type="checkbox" <?php if (isset($_GET['scale_atmosphere']) && $_GET['scale_atmosphere'] == 1) {?> checked="checked" <?php } ?> name="scale_atmosphere" value="1"> <b> <FONT COLOR="990000"> Check this box to scale the model by changing the scale height using actual temperature and gravity. Otherwise, the model will directly be used.</FONT></b>

		<?php
			if (isset($_GET['scale_atmosphere']) && $_GET['scale_atmosphere'] == 1) {
				$scale_atmosphere = 1;
			} else {
				$scale_atmosphere = 0;
			}
		?>
		</td></tr>
		<tr><td>Star Radius</td>
		<td><input type="text" name="star_radius" value="<?php echo isset($star_radius) ? $star_radius :  $_GET['star_radius'] ?>"  size="6" maxlength="5"/> Rsun</td></tr>
		<tr><td>Star Teff</td>
		<td><input type="text" name="star_teff" value="<?php echo isset($star_teff) ? $star_teff : $_GET['star_teff'] ?>"  size="6" maxlength="5"/> Kelvin (rounded to 100 K)</td></tr>
		<tr><td>Star J Mag</td>
		<td><input type="text" name="star_jvega" value="<?php echo isset($star_jvega) ? $star_jvega : $_GET['star_jvega'] ?>"  size="6" maxlength="5"/> (Vega System)</td></tr>
		<tr><td>Transit Duration</td>
		<td><input type="text" name="transit_duration" value="<?php echo isset($transit_duration) ? $transit_duration :  $_GET['transit_duration'] ?>"  size="6" maxlength="5"/> hours (t_14)</td></tr>
		<tr><td>Planet Density</td>
		<td><input type="text" name="planet_density" value="<?php echo isset($planet_density) ? $planet_density : $_GET['planet_density'] ?>"  size="6" maxlength="5"/> g/cm3</td></tr>
		<tr><td>Planet Teq</td>
		<td><input type="text" name="planet_teff" value="<?php echo isset($planet_teff) ? $planet_teff : $_GET['planet_teff'] ?>"  size="6" maxlength="4"/> Kelvin</td></tr>
		<tr><td>Planet Solid Radius</td>
		<td><input type="text" name="solid_planet_radius" value="<?php echo isset($solid_planet_radius) ? $solid_planet_radius :  $_GET['solid_planet_radius'] ?>"  size="6" maxlength="5"/>

		<?php
			$spr_unit_selected = isset($spr_unit) ? $spr_unit : $_GET['spr_unit'];
		?>
		<select name="spr_unit">
			<!--<option value="Rjup">Rjup</option>
			<option value="Rearth">Rearth</option>-->
  			<option <?php if($spr_unit_selected == "Rjup"){echo("selected");}?>>Rjup</option>
  			<option <?php if($spr_unit_selected == "Rearth"){echo("selected");}?>>Rearth</option>
		</select></td></tr>

		<tr><td colspan=2><b>Instrument Setup and Observing Strategy</b></td></tr>
		<tr><td>Out to In Factor</td>
		<td><input type="text" name="out_to_in_factor" value="<?php echo isset($_GET['out_to_in_factor']) ? $_GET['out_to_in_factor'] : $out_to_in_factor ?>"  size="6" maxlength="4"/> (Ratio of time spent outside transit and during transit)</td></tr>
		<tr><td>Number of Transits</td>
		<td><input type="text" name="ntransit" value="<?php echo isset($_GET['ntransit']) ? $_GET['ntransit'] : $ntransit ?>"  size="3" maxlength="3"/></td></tr>
		<tr><td>Pixel Binning  (Order 1)</td>
		<td><input type="text" name="pixelbin_m1" value="<?php echo isset($_GET['pixelbin_m1']) ? $_GET['pixelbin_m1'] : $pixelbin_m1 ?>"  size="3" maxlength="3"/></td></tr>
		<tr><td>Pixel Binning (Order 2)</td>
		<td><input type="text" name="pixelbin_m2" value="<?php echo isset($_GET['pixelbin_m2']) ? $_GET['pixelbin_m2'] : $pixelbin_m2 ?>"  size="3" maxlength="3"/> (2 for Nyquist sampling)</td></tr>
		<tr><td>Noise Floor (ppm)</td>
		<td><input type="text" name="floor_ppm" value="<?php echo isset($_GET['floor_ppm']) ? $_GET['floor_ppm'] : '30' ?>"  size="3" maxlength="3"/></td></tr>
		<tr><td>Sub-Array</td>

		<?php
			$readoutmode_selected = $_GET['readoutmode'];
		?>
		<td><select name="readoutmode">
			<!--<option value="256">256</option>
			<option value="96">96</option>
			<option value="2048">2048</option>-->
			<!--<option value="0">0</option>-->
 			<option <?php if($readoutmode_selected == "256"){echo("selected");}?>>256</option>
  			<option <?php if($readoutmode_selected == "96"){echo("selected");}?>>96</option>
  			<option <?php if($readoutmode_selected == "2048"){echo("selected");}?>>2048</option>
		</select>x2048</td></tr>
		<tr><td>Number of Read Frames per Ramp (aka Ngroup)</td>
		<td><input type="text" name="nread" value="<?php echo isset($_GET['nread']) ? $_GET['nread'] : $nread ?>"  size="3" maxlength="3"/></td></tr>
		<tr><td colspan=2><b><input type=submit value="Press to Run the Simulator" name="RunSimuButton"></b></td></tr>
	</table>
</form>
</div>

<?php $targetname = preg_replace(" ","", $targetname); ?>


<div id="p">
<!-- Ternary operator in PHP: echo expression ? $foo : $bar; 
Evaluates to the value of the second expression if the first one evaluates to TRUE, 
and evaluates to the third expression if the first evaluates to FALSE. -->

	<?php
	if (isset($_GET['RunSimuButton'])) {

		// Generate a unique string in case multiple queries simultaneously
		$bytes = rand(10000,99999);
		$runID_short = "runID_".$bytes;
		$runID = getcwd()."/public/".$runID_short."/";
		//echo $runID;
		//echo $runID_short;
		mkdir($runID,0777);
		chmod($runID,0777);

		if(isset($_GET['formPlanetModel']) ) {
			$planet_model = $_GET['formPlanetModel'];
		} else {
			$planet_model = "Not a valid planet model entered";
		}

		$myFile = $runID."simu1D_input.config";
		//echo $myFile;
		$fh = fopen($myFile, 'w') or die("Error. Can not read the simulation input config file. Send Loic a note...");
		$targetname = preg_replace('/\s+/','',$_GET['targetname']); //remove white spaces
		//echo $targetname;
		fwrite($fh, $targetname." ");
		fwrite($fh, $planet_model." ");
		fwrite($fh, $_GET['star_radius']." ");
		fwrite($fh, $_GET['star_teff']." ");
 		fwrite($fh, $_GET['star_jvega']." ");
 		fwrite($fh, $_GET['transit_duration']." "); 
		fwrite($fh, $_GET['planet_density']." "); 
		fwrite($fh, $_GET['planet_teff']." "); 
		fwrite($fh, $_GET['solid_planet_radius']." "); 
		fwrite($fh, $_GET['spr_unit']." "); 
		fwrite($fh, $scale_atmosphere." "); 
		fwrite($fh, $_GET['out_to_in_factor']." "); 
		fwrite($fh, $_GET['ntransit']." "); 
		fwrite($fh, $_GET['pixelbin_m1']." "); 
		fwrite($fh, $_GET['pixelbin_m2']." "); 
		fwrite($fh, $_GET['floor_ppm']." "); 
		fwrite($fh, $_GET['readoutmode']." "); 
		if($_GET['nread'] == 0) {
		// carefull here - remove this once I implement the algorithm that sets the NREAD automatically.
			$nread = 2;
		} else {
			$nread = $_GET['nread'];
		}
		fwrite($fh, $nread." ");
		fclose($fh);
	} else {
		echo "Make sure to fill in the boxes before pressing the button.Running the simulation takes about 2 seconds. Patience. ";
	}

	// Create a ghost file for testing
	$myFile_new = $runID."simu1D_input.config";
	$fh_new = fopen($myFile_new, 'w') or die("...");
	$c1 = $targetname;
	$c2 = $planet_model;
	$c3 = $_GET['star_radius'];
	$c4 = $_GET['star_teff'];
	$c5 = $_GET['star_jvega'];
	$c6 = $_GET['transit_duration'];
	$c7 = $_GET['planet_density'];
	$c8 = $_GET['planet_teff'];
	$c9 = $_GET['solid_planet_radius'];
	$c10 = $_GET['spr_unit'];
	$c11 = $scale_atmosphere;
	$c12 = $_GET['out_to_in_factor'];
	$c13 = $_GET['ntransit'];
	$c14 = $_GET['pixelbin_m1'];
	$c15 = $_GET['pixelbin_m2'];
	$c16 = $_GET['floor_ppm'];
	$c17 = $_GET['readoutmode'];
	$c18 = $nread;
	$line = sprintf("%-31s %-60s %-6s %6s %-6s %-6s %-6s %6s %-6s %6s %3s %5s %4s %4s %4s %4s %5s %4s\n",$c1,$c2,$c3,$c4,$c5,$c6,$c7,$c8,$c9,$c10,$c11,$c12,$c13,$c14,$c15,$c16,$c17,$c18);
	fwrite($fh_new, $line);
	fclose($fh_new);

	//if (is_file($myfile)) {
		//worked before $command = "/usr/local/bin/idl -e \"simu1d,"." configfile='".$runID."simu1D_input.config'\""." > ".$runID."out.txt 2> ".$runID."err.txt";
		//calling pref  $command = "/usr/local/bin/idl -pref=\"/home/cpapir/www/niriss/simu1D/idlpath.txt\" -e \"simu1d,"." configfile='".$runID."simu1D_input.config'\""." > ".$runID."out.txt 2> ".$runID."err.txt";
      	exec("/usr/local/bin/idl runsimu1d.pro -args '".$runID."simu1D_input.config'"." > ".$runID."out.txt 2> ".$runID."err.txt");

		echo $command;
		exec($command);

	//}
	//change mode of created files
	chmod($runID."simu1D_input.config",0666);
	chmod($runID."out.txt",0666);
	chmod($runID."err.txt",0666);
	chmod($runID."out2.txt",0666);
	chmod($runID."err2.txt",0666);
	chmod($runID."simu1D_results_order1.txt",0666); 
	chmod($runID."simu1D_results_order2.txt",0666); 
	chmod($runID."simu1D_parameters.txt",0666);
	?>

	<?php
	//$makefig_command = getcwd()."/simu1D_fig.py ".$targetname." ".$runID."simu1D_results_order1.txt ".$runID."simu1D_results_order2.txt ".$runID."simu1D_parameters.txt "." > ".$runID."out2.txt 2> ".$runID."err2.txt";
	$makefig_command = getcwd()."/simu1D_fig.py ".$targetname." ".$runID."simu1D_results_order1.txt ".$runID."simu1D_results_order2.txt ".$runID."simu1D_model_order1.txt ".$runID."simu1D_model_order2.txt ".$runID."simu1D_parameters.txt "." > ".$runID."out2.txt 2> ".$runID."err2.txt";

	//echo $makefig_command;
	exec($makefig_command);
	?>

	<?php
	$maketarball = "tar -cvf ".$runID.$runID_short."_tarball.tar -C ".getcwd()."/public/".$runID_short." .";
	//echo $maketarball;
	exec($maketarball);
	$gziptarball = "gzip ".$runID.$runID_short."_tarball.tar";
	//echo $gziptarball;
	exec($gziptarball);
	?>

	<div id="p">
	<h2>Table of Simulation Results</h2>
	<table align=left border="1">
	<tr>
	<td><b>Input Parameters:</a></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/simu1D_parameters.txt">simu1D_parameters.txt</a></b></td>
	<td></td>
	</tr>
	<tr>
	<td><b>ASCII Output:</a></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/simu1D_results_order1.txt">simu1D_results_order1.txt</a></b></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/simu1D_results_order2.txt">simu1D_results_order2.txt</a></b></td>
	</tr>
	<tr>
	<td></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/simu1D_model_order1.txt">simu1D_model_order1.txt</a></b></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/simu1D_model_order2.txt">simu1D_model_order2.txt</a></b></td>
	</tr>
	<tr>
	<td><b>Figures:</a></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/PhotonCount_<?php echo $targetname ?>.pdf">Photon Count Plot (pdf)</a></b></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/Signal2Noise_<?php echo $targetname ?>.pdf">Signal to Noise Plot (pdf)</a></b></td>
	</tr>
	<tr>
	<td></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/NoisePPM_<?php echo $targetname ?>.pdf">Noise Plot (pdf)</a></b></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/ResolvingPower_<?php echo $targetname ?>.pdf">Resolving Power (pdf)</a></b></td>
	</tr>
	<tr>
	<td></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/SpectrumRealizationPPM_<?php echo $targetname ?>.pdf">Spectrum Realization in ppm (pdf)</a></b></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/TransitDepthPPM_<?php echo $targetname ?>.pdf">Transit Depth in ppm (pdf)</a></b></td>
	</tr>
	<tr>
	<td><b>Diagnostics:</a></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/out.txt">log output (calcs)</a></b></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/err.txt">log errors (calcs)</a></b></td>
	</tr>
	<tr>
	<td></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/out2.txt">log output (plots)</a></b></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/err2.txt">log errors (plots)</a></b></td>
	</tr>
	<tr>
	<td><b>Tar Ball of all Files:</a></td>
	<td><b><a href="<?php echo $url?>public/<?php echo $runID_short ?>/<?php echo $runID_short ?>_tarball.tar.gz">TarBall</a></b></td>
	<td></td>
	</tr>
	</table>



</div>


</body>

</html>
