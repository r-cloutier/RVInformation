<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<link rel="stylesheet" type="text/css" href="style.css" />
<title>Transiting RV Calculator</title>

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
        	<h1>Transiting RV Calculator</h1>
            <h2>Ryan Cloutier</h2>
        </div>   
        
        <div id="menu">
        	<!--<ul>
				<li class="menuitem"><a href="http://jgagneastro.wordpress.com/">Research Material</a></li>
            	<li class="menuitem"><a href="http://www.exoplanetes.umontreal.ca/?lang=en">Research Group</a></li>
                <li class="menuitem"><a href="http://www.astro.umontreal.ca/~malo/banyan.php">BANYAN I <br> web tool</a></li>
				<li class="menuitem"><a href="http://jgagneastro.wordpress.com/banyanii/">BANYAN II <br>Additional Material</a></li>
				<li class="menuitem"><a href="http://www.das.uchile.cl/~drodrigu/CPCalc.html">Convergent Point (David Rodriguez)</a></li>
            </ul> -->
        </div>
        

		<div id="content">
		  <h2>Calculate some stuff</h2>
		  
		  <div id="abstract">
		    <br><b> Talk about what the damn tool does.</b><br>
		  </div>

		  <form action="TransitRVCalculator.php" method="get" >
		    <a1><br>Param1 : <input type="text" name="param1" value="<?php echo isset($_GET['param1']) ? $_GET['param1'] : $param1 ?>" size="10" maxlength="50"/><br>
		      <br><input type=submit value="Submit dat" name="submit"/><br>
		      <br>  
		      <br></a1>
		  </form>

		</div>


		<div>
		  <?php
		     if (isset($_GET['submit'])) {

		     echo $_GET["p1"];
		     
		     }

		     ?>

		</div>
		  
</body>

</html>
