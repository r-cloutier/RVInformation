<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<link rel="stylesheet" type="text/css" href="style.css" />
<title>BANYAN II </title>

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

</head>

<body>

<div id="container">
		<div id="mainpic">
        	<h1>BANYAN II web tool</h1>
            <h2>by Jonathan Gagné</h2>
        </div>   
        
        <div id="menu">
        	<ul>
				<li class="menuitem"><a href="http://jgagneastro.wordpress.com/">Research Material</a></li>
            	<li class="menuitem"><a href="http://www.exoplanetes.umontreal.ca/?lang=en">Research Group</a></li>
                <li class="menuitem"><a href="http://www.astro.umontreal.ca/~malo/banyan.php">BANYAN I <br> web tool</a></li>
				<li class="menuitem"><a href="http://jgagneastro.wordpress.com/banyanii/">BANYAN II <br>Additional Material</a></li>
				<li class="menuitem"><a href="http://www.das.uchile.cl/~drodrigu/CPCalc.html">Convergent Point (David Rodriguez)</a></li>
<!--				<li class="menuitem"><a href="http://figshare.com/authors/Jonathan%20Gagn%C3%A9/479804">J. Gagné's<br> Figshare Material</a></li> -->
<!--				<li class="menuitem"><a href="https://www.researchgate.net/profile/Jonathan_Gagne/contributions/?ev=prf_act">J. Gagné's<br> ResearchGate Material</a></li> -->
            </ul>
        </div>
        

		<div id="content">
			<h2>BANYAN II: Bayesian Analysis for Nearby Young AssociatioNs II</h2> 
			<h1>Membership probability without photometric information </h1> 
			<h1>(Version 1.4)</h1> 
		
			<div id="abstract">
				<br><b> BANYAN II is a Bayesian analysis tool to determine the membership probability of candidate 
					stars to nearby young kinematic groups. This pipeline is based on a comparison of Galactic 
					position (XYZ) and space velocity (UVW), using a naive Bayesian classifier. We only consider 
					well-defined moving groups closer than 100 pc and younger than 200 Myr. When radial velocity
					and/or distance measurements are not available, Bayesian inference has the advantage of being
					able to marginalize them and compute membership probabilities without them. When this is the case,
					BANYAN II produces statistical predictions for these quantities, according to each membership
					hypothesis. When a given hypothesis is respected, the agreement between trigonometric distances
					and BANYAN II statistical distances generally agree within 8%, whereas predicted and measured radial
					velocities agree within 1.6 km/s. <br><br>
					
					All modifications included in this BANYAN II web tool with respect to the 
					<a href="http://www.astro.umontreal.ca/~malo/banyan.php">BANYAN I web tool</a> are listed here :<br><br>
					
					<li> The spatial and kinematic models of moving groups are modeled with 3D gaussian ellipsoids in both
						BANYAN I and II, but in BANYAN II those ellipsoids are free to have axes oriented in any direction,
						rather than being constrained in the local Galactic coordinates (XYZ) directions. </li><br>
				    <li> The list of bona fide members used in the construction of the spatial and kinematic models
					    was updated as of December 2013.</li><br>
					<li> The field hypothesis is modeled using a Besançon Galactic model (Robin et al. 2012).</li><br>				
					<li> Prior probabilities are not set to unity anymore, but rather to the expected populations in
						each hypothesis considered, to bring the Bayesian probabilities closer to absolute values.
						This treatment of relative populations takes into account that objects with A) smaller galactic
						latitude, b) smaller proper motion, c) smaller RVs and d) larger distances, are more likely to
						be field contaminants. </li><br>					
					<li> The likelihood functions that are input in Baye's theorem are functions of X,Y,Z,U,V,W instead of
						the direct observable quantities (ra, dec, pmra, pmdec, RV, distance). The reason for this preference
						is that former quantities are better represented by gaussian likelihood probability density functions. </li><br>
					<li> Prior likelihoods associated to marginalized parameters (distance and/or RV) are computed directly
						from a Montecarlo drawing of synthetic objects from each association's spatial and kinematic model,
						which means that they can have any functional form (generally not gaussians). </li><br>
					<br>
					The details of this method are described in <a href="http://adsabs.harvard.edu/abs/2014ApJ...783..121G">Gagné et al. (2014ApJ...783..121G)</a>. Please reference
					this work when you make use of this tool, as well as <a href="http://adsabs.harvard.edu/abs/2013ApJ...762...88M"> Malo et al. (2013)</a>. </b><br><br>
					
					We also encourage trying <a href="http://www.das.uchile.cl/~drodrigu/CPCalc.html">David Rodriguez' convergent point analysis online tool</a>.<br><br>

					We would like to thank Michael Liu for useful comments on this web page.<br>

					<br> <b><FONT SIZE="3">Version history</FONT></b>.<br><br>

						<b>v1.0</b> <i>22/12/2013</i><br>
						<b>v1.1</b> <i>24/12/2013</i> Probabilities taking into account only a RV or PLX measurement (but not both) were mistaken, due to a misdistribution of prior probabilities (every hypothesis always inherited the Beta Pictoris prior probability - depending on some cases, this could have a very small or very big effect on the probabilities). This has been fixed from this version, as well as in the v2+ of the arxiv version of Gagné et al. (2014). The ApJ version was unaffected.<br>
						<b>v1.2</b> <i>3/1/2014</i> Probabilities taking into account either a RV or PLX measurement (or both) were slightly offset, because the &xi;<sub>&nu;</sub> (related to RV) and &xi;<sub>&varpi;</sub> (related to distance) were always = 1 (See Section 3.1 of Gagné et al. 2014 for more info). This has been fixed from this version. This error was not present in the paper, and thus did not affect either the arxiv or ApJ version.<br>
						<b>v1.3</b> <i>23/4/2014</i> RV and PLX entries without decimals were previously (mistakenly) treated as integers, hence yielding problematic results.<br>
						<b>v1.4</b> <i>11/8/2015</i> Separated the output probabilities of "Young" and "Old" field, added the statistical RVs for the field hypotheses, and fixed glitches in the output format of statistical RVs and distances.<br>


				<br> <b><FONT COLOR="FF0000" SIZE="5">The results of this tool should be interpreted with <u>caution</u></FONT></b>.<br><br>
				
					<b>1)</b> We stress that a high 
					probability in one of the associations considered here does not necessarily imply that the 
					candidate star is young. Before an object is considered as a new bona fide member to a moving group, 
					it must show signs of youth as well as spatial and kinematic agreement with the moving group in question
					(using both a RV and parallax measurements). Conversely, a high probability in the &laquo; field &raquo; 
					hypothesis does not necessarily imply that the star is old; it simply means that it is less likely to be 
					a member of the young associations considered here compared to the &laquo; field &raquo;<br><br>
					
					<b>2)</b> Users must keep in mind that Bayesian probabilities might be wrong when studying members of
					moving groups or associations not considered in our hypotheses. For example, a member of Upper Scorpius
					could come up as a good candidate of AB Doradus, because the models used by our tool do not consider the
					former.<br><br>

					<b>3)</b> Bayesian probabilities reported by BANYAN II will often be lower than those reported by BANYAN I,
					because the priors were not set to unity here, which strongly favors the field that generally has a significantly
					higher expected population. This means in turn that any object with Bayesian probabilities not negligible for a given
					moving group is worthy of further study (e.g. a 20% probability in BANYAN II could amount to a 90% probability in BANYAN I).<br><br>

					<b>4)</b> Bayesian probabilities calculated using a naive Bayesian classifier (which is the case in BANYAN I and II)
					are subject to being biased because of the independent treatment of observables that are not independent in reality.
					For a given object, this will not affect the relative classification rank of hypotheses (e.g. from least probable to
					most probable), but it will affect the absolute values of probabilities (see Hand and Yu 2001).<br><br>

					<b>5)</b> Results obtained here may differ from those published in <a href="http://adsabs.harvard.edu/abs/2014ApJ...783..121G">Gagné et al. (2014ApJ...783..121G)</a> as the latter includes 2MASS
					and WISE photometry in computing membership probabilities.<br><br>
					
					<b>6)</b> The information entered on this web page is not stored on our server. <br><br>
			
			<br> <b><FONT SIZE="5">Instructions :</FONT></b><br><br>          
			
	    		<li> You can simply enter the name of a star, as resolved by Simbad or Vizier.</li><br>
	    		<li> Press the RESOLVE button. BANYAN II will seek all the information it can find from a limited number 
					of online catalogs. Missing information will be reported as &laquo; NaN &raquo;.Be patient, this may 
					take several tens of seconds.</li><br>
	    		<li> You can modify and/or remove some information as you please but you need minimally a sky position, proper motion
					and a measurement error on proper motion to proceed.</li><br>
				<li> Proper motions are in units of mas/yr, and the proper motion in the RA direction is implicitly understood as the usual mu_ra * cos(delta), where delta is the declination of the star.</li><br>
          		<li> Press SUBMIT and read the cautionary note above (if not already done) before interpreting the results.</li><br>
          		<li> You can use the "set priors to unity" checkbox in order to compare probabilities with BANYAN I. However, these probabilities will be more strongly biased (far from absolute probabilities ; See Gagné et al. 2014.)</li><br>
          		<li> The main effect of checking "Object is younger than 1 Gyr" is that prior probabilities are different (apart from the fact that the field population dymanics is dependent on age) - this is so because the expected population of old field objects is significantly larger than that of young objects. Hence, this option will generally have a very small effect on Bayesian probabilities if the "set priors to unity" option is checked.</li><br>
				<li> If you check the "Younger than 1 Gyr" option, the "Old field" hypothesis is withdrawn from the analysis, and hence its associated probability is automatically zero. In that case, all the reported probabilities, statistical distances and statistical RVs are for the young field. If you do not check the box, separate probabilities, RV and distance predictions will be given out - "OFLD" will correspond to the old field, and "YFLD" to the young field.</li><br>
				<li> Please be reminded that this tool does not use photometry as input observables. Most of the papers about brown dwarfs in the BANYAN series use an augmented version of BANYAN II that does use 2MASS and WISE photometry, hence the resulting output probabilities will be different.</li><br>
          		<li> Acknowledgement : if your paper uses results obtained with BANYAN II, please kindly give a reference to <a href="http://adsabs.harvard.edu/abs/2014ApJ...783..121G">Gagné et al. (2014ApJ...783..121G)</a> and <a href="http://adsabs.harvard.edu/abs/2013ApJ...762...88M"> Malo et al. (2013)</a>. </li><br>
			    <li> If you have any questions and/or comments, contact me at gagne (@) astro (dot) umontreal (dot) ca.</li><br>
			
			</div>

			<form action="http://www.astro.umontreal.ca/~gagne/banyanII.php" method="get" >
			NAME of your star: <input type="text" name="targetname" value="<?php echo isset($_GET['targetname']) ? $_GET['targetname'] : '' ?>"  size="30" maxlength="110"/>
   			PRESS: <input type=submit value="Resolve" name="resolve"><br>
			</form>

			<?php
				if (isset($_GET['resolve'])) {

				$dd = $_GET['targetname'];
				$name = "\"&".$dd."\"";
				$name = $name;

				$row = 0;
				$pp = "/home/gagne/www/banya/answer/info_";
				$pp1 = ".dat";
				$pp2 = $pp.$dd.$pp1;
				$Fichier = $pp2;

				if (is_file($Fichier)) {
				} else {
				$command = "/usr/local/itt/idl/idl/bin/idl -e 'jonathan_name_banyanii,name=$name' ";
				exec($command);
				}
	
				if (($handle = fopen($Fichier, "r")) !== FALSE) {
					while (($data = fgetcsv($handle, 1000, ",")) !== FALSE) {
    					if ($row == 0) { 
       				 // this is the first line of the csv file
        				// it usually contains titles of columns
        				$num = count($data);
	  				$row++;
    					} else {
						$ra = $data[1];
						$dec = $data[2];
						$pmra = $data[3];
						$epmra = $data[4];
						$pmdec = $data[5];
						$epmdec = $data[6];
						$refpm = $data[7];
						$vrad = $data[8];
						$evrad = $data[9];
						$plx = $data[10];
						$eplx = $data[11];
   						 }
					}
					fclose($handle);
				}

				}
				
			?>
				
				<form action="http://www.astro.umontreal.ca/~gagne/banyanII.php" method="get" >
				<a1><br><input type="checkbox" <?php if (isset($_GET['isyoung']) && $_GET['isyoung'] == 'Yes'){?> checked="checked" <?php } ?> name="isyoung" value="Yes"> <b> <FONT COLOR="990000"> Check this box if this object is younger than 1 Gyr.</FONT></b>
				<br><input type="checkbox" <?php if (isset($_GET['noprior']) && $_GET['noprior'] == 'Yes'){?> checked="checked" <?php } ?> name="noprior" value="Yes"> <b> <FONT COLOR="990000"> Check this box to set prior probabilities to unity.</FONT></b>
				<br><br>Right ascension (degrees) : <input type="text" name="radeg" value="<?php echo isset($_GET['radeg']) ? $_GET['radeg'] : $ra ?>"  size="10" maxlength="50"/> <br>
		 		<br>Proper motion in right ascension (mas/yr) : <input type="text" name="pmra" value="<?php echo isset($_GET['pmra']) ? $_GET['pmra'] : $pmra ?>"  size="10" maxlength="50"/><br>
 				 <br>Proper motion in declination (mas/yr) : <input type="text" name="pmdec" value="<?php echo isset($_GET['pmdec']) ? $_GET['pmdec'] : $pmdec ?>"  size="10" maxlength="50"/><br>
				 <br>Radial velocity (km/s) : <input type="text" name="hrv" value="<?php echo isset($_GET['hrv']) ? $_GET['hrv'] : $vrad ?>"  size="10" maxlength="50"/><br>
				 <br>Parallax (mas) : <input type="text" name="plx" value="<?php echo isset($_GET['plx']) ? $_GET['plx'] : $plx  ?>"  size="10" maxlength="50"/><br>		
				 <br> STEP 2 : <input type=submit value="Submit" name="submit"/><br>
				 <br>  
				 <br>  </a1>

 				<a2><br><br><br><br>Declination (degrees) : <input type="text" name="dec" value="<?php echo isset($_GET['dec']) ? $_GET['dec'] : $dec ?>"  size="10" maxlength="50"/><br>
				 <br>Error on Proper motion in right ascension  (mas/yr) : <input type="text" name="epmra" value="<?php echo isset($_GET['epmra']) ? $_GET['epmra'] : $epmra ?>"  size="10" maxlength="50"/><br>
 				 <br>Error on Proper motion in declination  (mas/yr) : <input type="text" name="epmdec" value="<?php echo isset($_GET['epmdec']) ? $_GET['epmdec'] : $epmdec ?>"  size="10" maxlength="50"/><br>
				 <br>Error on radial velocity (km/s) : <input type="text" name="ehrv" value="<?php echo isset($_GET['ehrv']) ? $_GET['ehrv'] : $evrad  ?>"  size="10" maxlength="50"/><br>
				 <br>Error on parallax (mas) : <input type="text" name="eplx" value="<?php echo isset($_GET['eplx']) ? $_GET['eplx'] : $eplx  ?>"  size="10" maxlength="50"/><br>		
				 <br> Name of your star: <input type="text" name="targetname" value="<?php echo isset($_GET['targetname']) ? $_GET['targetname'] : $name ?>"  size="30" maxlength="50"/><br>
				 <br>  </a2>
				</form>

			<?php

				if (isset($_GET['submit'])) {
				
				$dd = $_GET['targetname'];
				$name = "\"&".$dd."\"";
				$name2 = "\"".$dd."\"";
				$ra = $_GET['radeg'];
				$dec = $_GET['dec'];
				$pmra = $_GET['pmra'];
				$epmra = $_GET['epmra'];
				$pmdec = $_GET['pmdec'];
				$epmdec = $_GET['epmdec'];
				$vrad = $_GET['hrv'];
				$evrad = $_GET['ehrv'];
				$plx = $_GET['plx'];
				$eplx = $_GET['eplx'];
				
				if(isset($_GET['isyoung']) && $_GET['isyoung'] == 'Yes')
				{
				    $isyoung = "0";
				}
				else
				{
				    $isyoung = "1";
				}
				if(isset($_GET['noprior']) && $_GET['noprior'] == 'Yes')
				{
				    $noprior = "1";
				}
				else
				{
				    $noprior = "0";
				}
				
					if (!is_numeric($vrad))
					{
					$vrad=-99.00;
					} 
					if (!is_numeric($evrad))
					{
					$evrad=-99.00;
					} 

					if (!is_numeric($plx))
					{
					$plx=-99.00;
					} 

					if (!is_numeric($eplx))
					{
					$eplx=-99.00;
					} 
				$ip = $_SERVER['REMOTE_ADDR'];
				$d = "\n";
				$ips = "\"&".$ip."\"";
				$command = "/usr/local/itt/idl/idl/bin/idl -e 'jonathan_asso_banyanii,name=$name,ra=$ra,dec=$dec,pmra=$pmra,epmra=$epmra,pmdec=$pmdec,epmdec=$epmdec,vrad=$vrad,evrad=$evrad,plx=$plx,eplx=$eplx,old=$isyoung,noprior=$noprior,ip=$ips' ";
				exec($command);
				
				echo "<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br> The membership probabilities (%) for $name2 are : <br>";

			$myfile = "banya/answer/submit_visitors.txt";
			$fh = fopen($myfile, 'a');
			$stringData = $ip.$d;
			fwrite($fh, $stringData);
			fclose($fh);

			$count_my_page = "banya/answer/submit_hit.txt";
			$hits = file($count_my_page);
			$hits[0] ++;
			$fp = fopen($count_my_page , "w");
			fputs($fp , "$hits[0]");
			fclose($fp);

			?>

			<table border="1">
				<?php

					$row = 0;
					$pp = "/home/gagne/www/banya/answer/prob_";
					$pp1 = ".dat";
					$pp2 = $pp.$dd.$pp1;
					$Fichier = $pp2;

					if (($handle = fopen($Fichier, "r")) !== FALSE) {
						while (($data = fgetcsv($handle, 1000, ",")) !== FALSE) {
    							if ($row == 0) { 
        						// this is the first line of the csv file
       						 // it usually contains titles of columns
       						 $num = count($data);
	
        						echo "<thead>\n<tr style=\"background-color:#CC0000;\">";
        						$row++;
        						for ($c=0; $c < $num; $c++) { 
            						echo "<th><font size='2'>" . $data[$c] . "</font></th>";
       						 }
       							 echo "</tr>\n</thead>\n\n<tbody>";
    							} else {
        						// this handles the rest of the lines of the csv file
        						$num = count($data);
        						echo "<tr>";
        						$row++;
        						for ($c=0; $c < $num; $c++) {
           						 echo "<td><font size='2'>" . $data[$c] . "</font></td>";
        						}
        						echo "</tr>\n";
    							}
						}
					fclose($handle);
				}

				?>
			</table>

			<?php
				echo "<b>Legend</b>: P_ASSOCIATION: probability which takes into account RA,DEC and proper motion. 
				<br> PV_ASSOCIATION: probability which takes into account RA,DEC, proper motion and radial velocity.
				<br> PP_ASSOCIATION: probability which takes into account RA,DEC, proper motion and parallax.
				<br> PVP_ASSOCIATION:  probability which takes into account RA,DEC, proper motion, radial velocity and parallax.<br> " ;
			?>


			<table border="1">
				<?php
				echo "<BR> The predicted radial velocities for $name2 are :";

				$row = 0;
				$pp = "/home/gagne/www/banya/answer/mvrad_";
				$pp1 = ".dat";
				$pp2 = $pp.$dd.$pp1;
				$Fichier = $pp2;

				if (($handle = fopen($Fichier, "r")) !== FALSE) {
					while (($data = fgetcsv($handle, 1000, ",")) !== FALSE) {
    						if ($row == 0) { 
        					// this is the first line of the csv file
        					// it usually contains titles of columns
        					$num = count($data);
	
        				echo "<thead>\n<tr style=\"background-color:#CC0000;\">";
        				$row++;
        					for ($c=0; $c < $num; $c++) { 
            				echo "<th><font size='1'>" . $data[$c] . "</font></th>";
        					}
        					echo "</tr>\n</thead>\n\n<tbody>";
    						} else {
        					// this handles the rest of the lines of the csv file
        					$num = count($data);
        					echo "<tr>";
        					$row++;
        						for ($c=0; $c < $num; $c++) {
            					echo "<td><font size='2'>" . $data[$c] . "</font></td>";
       						 }
        						echo "</tr>\n";
    							}
						}
					fclose($handle);
				}

				?>
			</table>


			<table border="1">
			<?php
				echo "<BR> The predicted distances for $name2 are :";

				$row = 0;
				$pp = "/home/gagne/www/banya/answer/mdist_";
				$pp1 = ".dat";
				$pp2 = $pp.$dd.$pp1;
				$Fichier = $pp2;

				if (($handle = fopen($Fichier, "r")) !== FALSE) {
					while (($data = fgetcsv($handle, 1000, ",")) !== FALSE) {
    						if ($row == 0) { 
        					// this is the first line of the csv file
        					// it usually contains titles of columns
        					$num = count($data);
	
        					echo "<thead>\n<tr style=\"background-color:#CC0000;\">";
        					$row++;
        						for ($c=0; $c < $num; $c++) { 
            					echo "<th><font size='1'>" . $data[$c] . "</font></th>";
        						}
        						echo "</tr>\n</thead>\n\n<tbody>";
    							} else {
        						// this handles the rest of the lines of the csv file
        						$num = count($data);
        						echo "<tr>";
        						$row++;
        						for ($c=0; $c < $num; $c++) {
            					echo "<td><font size='2'>" . $data[$c] . "</font></td>";
        						}
        						echo "</tr>\n";
    							}
						}
				fclose($handle);
				}
				
			?>
			</table>


			<?php
			}
				echo "<center> <br><br><br>Observational properties (UVWXYZ) of the young associations considered in the analysis : </center>";
				echo "<center>     *Note that primed quantities are in measured in a rotated frame of reference (see Gagne et al. 2014) : </center>";
			?>

			<div id="donnee">

<table border="1">
<th> Name</th>
<th> Age (Myr)</th>
<th> 1-&sigma; distance range (pc)</th>
<th> 1-&sigma; RV range (km/s)</th>
<th> U (km/s)</th>
<th> V (km/s)</th>
<th> W (km/s)</th>
<th> &sigma; U' (km/s)</th>
<th> &sigma; V' (km/s)</th>
<th> &sigma; W' (km/s)</th>
<th> X (pc)</th>
<th> Y (pc)</th>
<th> Z (pc)</th>
<th> &sigma; X' (pc)</th>
<th> &sigma; Y' (pc)</th>
<th> &sigma; Z' (pc)</th>
<th> # stars </th>

<tr>
	<td>TW Hydrae</td>
	<td> 8 - 12 </td>
	<td> 40 - 62 </td>
	<td> (7) - (12) </td>
	<td> -11.1 </td>
	<td> -18.9 </td>
	<td> -5.6 </td>
	<td> 0.9 </td>
	<td> 1.6 </td>
	<td> 2.8 </td>
	<td> 19.1 </td>
	<td> -54.2 </td>
	<td> 21.5 </td>
	<td> 5.0 </td>
	<td> 7.2 </td>
	<td> 22.6 </td>
	<td> 18 </td>	
</tr>
<tr>
	<td>&beta; Pictoris</td>
	<td> 12 - 22 </td>
	<td> 18 - 40 </td>
	<td> (-9) - (16) </td>
	<td> -11.0 </td>
	<td> -15.6 </td>
	<td> -9.2 </td>
	<td> 1.4 </td>
	<td> 1.7 </td>
	<td> 2.5 </td>
	<td> 7.6 </td>
	<td> -3.5 </td>
	<td> -14.5 </td>
	<td> 8.2 </td>
	<td> 13.5 </td>
	<td> 30.7 </td>
	<td> 33 </td>
</tr>
<tr>
<td>Tucana-Horologium</td>
	<td> 20 - 40 </td>
	<td> 38 - 51 </td>
	<td> (3) - (14) </td>
	<td> -9.7 </td>
	<td> -20.5 </td>
	<td> -0.8 </td>
	<td> 1.1 </td>
	<td> 1.7 </td>
	<td> 2.4 </td>
	<td> 6.7 </td>
	<td> -21.8 </td>
	<td> -36.0 </td>
	<td> 3.9 </td>
	<td> 10.6 </td>
	<td> 20.1 </td>
	<td> 52 </td>
</tr>
<tr>
	<td>Columba</td>
	<td> 20 - 40 </td>
	<td> 26 - 63 </td>
	<td> (19) - (26) </td>
	<td> -12.1 </td>
	<td> -21.3 </td>
	<td> -5.6 </td>
	<td> 0.5 </td>
	<td> 1.3 </td>
	<td> 1.7 </td>
	<td> -28.1 </td>
	<td> -25.8 </td>
	<td> -28.6 </td>
	<td> 10.5 </td>
	<td> 17.6 </td>
	<td> 28.3 </td>
	<td> 21 </td>
</tr>
<tr>
	<td>Carina</td>
	<td> 20 - 40 </td>
	<td> 11 - 42 </td>
	<td> (16) - (23) </td>
	<td> -10.7 </td>
	<td> -22.2 </td>
	<td> -5.7 </td>
	<td> 0.3 </td>
	<td> 0.7 </td>
	<td> 1.1 </td>
	<td> 10.1 </td>
	<td> -51.6 </td>
	<td> -14.9 </td>
	<td> 5.8 </td>
	<td> 11.3 </td>
	<td> 29.8 </td>
	<td> 21 </td>
</tr>
<tr>
	<td>Argus</td>
	<td> 30 - 50 </td>
	<td> 15 - 48 </td>
	<td> (-10) - (9) </td>
	<td> -21.5 </td>
	<td> -12.2 </td>
	<td> -4.6 </td>
	<td> 0.9 </td>
	<td> 1.7 </td>
	<td> 2.7 </td>
	<td> 15.0 </td>
	<td> -21.7 </td>
	<td> -8.1 </td>
	<td> 12.1 </td>
	<td> 15.5 </td>
	<td> 27.4 </td>
	<td> 11 </td>
</tr>
<tr>
	<td>AB Doradus</td>
	<td> 70 - 120 </td>
	<td> 19 - 50 </td>
	<td> (-11) - (29) </td>
	<td> -7.0 </td>
	<td> -27.2 </td>
	<td> -13.9 </td>
	<td> 1.2 </td>
	<td> 1.7 </td>
	<td> 1.9 </td>
	<td> -2.5 </td>
	<td> 1.3 </td>
	<td> -16.3 </td>
	<td> 16.3 </td>
	<td> 20.0 </td>
	<td> 23.5 </td>
	<td> 54 </td>
</tr>
<tr>
	<td>Young Field</td>
	<td> &#60; 1000 </td>
	<td> 66 - 169 </td>
	<td> (-19) - (19) </td>
	<td> -11.2 </td>
	<td> -18.6 </td>
	<td> -6.9 </td>
	<td> 7.7 </td>
	<td> 12.5 </td>
	<td> 19.6 </td>
	<td> 2.8 </td>
	<td> 0.1 </td>
	<td> -13.1 </td>
	<td> 79.6 </td>
	<td> 80.4 </td>
	<td> 80.8 </td>
	<td> - </td>
</tr>
<tr>
	<td>Old Field</td>
	<td> &#62; 1000 </td>
	<td> 70 - 177 </td>
	<td> (-34) - (32) </td>
	<td> -11.2 </td>
	<td> -18.6 </td>
	<td> -6.9 </td>
	<td> 7.7 </td>
	<td> 12.5 </td>
	<td> 19.6 </td>
	<td> 2.8 </td>
	<td> 0.1 </td>
	<td> -13.1 </td>
	<td> 79.6 </td>
	<td> 80.4 </td>
	<td> 80.8 </td>
	<td> - </td>
</tr>
</table>

</div>

<style>
table { text-align: center; border-collapse: collapse; font-size: 12px; }
tr:hover { background: blue; color: white }
th, td { padding: 5px }
</style>

            <div id="footer"><h3><a href="http://www.bryantsmith.com">florida web design</a></h3></div>


<?php
$count_my_page = ("/home/gagne/www/banya/answer/hitcounter_b2.txt");
$hits = file($count_my_page);
$hits[0] ++;
$fp = fopen($count_my_page , "w");
fputs($fp , "$hits[0]");
fclose($fp);
?>

<?php
$pp = "/home/gagne/www/banya/answer/prob_";
$ppp = "/home/gagne/www/banya/answer/mvrad_";
$pppp = "/home/gagne/www/banya/answer/mdist_";
$ppppp = "/home/gagne/www/banya/answer/info_";
$pp1 = ".dat";
$pp2 = $pp.$dd.$pp1;
$pp3 = $ppp.$dd.$pp1;
$pp4 = $pppp.$dd.$pp1;
$pp5 = $ppppp.$dd.$pp1;
unlink($pp2);
unlink($pp3);
unlink($pp4);
unlink($pp5);
$d1 = "/home/gagne/www/banya/answer/input.dat";
$d2 = "/home/gagne/www/banya/answer/output.dat";
unlink($d1);
unlink($d2);
?>

</div>
          
</div>

<!-- Start of StatCounter Code for Default Guide -->
<script type="text/javascript">
var sc_project=9506622; 
var sc_invisible=0; 
var sc_security="e913d2e1"; 
var scJsHost = (("https:" == document.location.protocol) ?
"https://secure." : "http://www.");
document.write("<sc"+"ript type='text/javascript' src='" +
scJsHost+
"statcounter.com/counter/counter.js'></"+"script>");
</script>
<noscript><div class="statcounter"><a title="web analytics"
href="http://statcounter.com/" target="_blank"><img
class="statcounter"
src="http://c.statcounter.com/9506622/0/e913d2e1/0/"
alt="web analytics"></a></div></noscript>
<!-- End of StatCounter Code for Default Guide -->

</body>

</html>
