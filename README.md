<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en"><head>
<meta name="generator" content="jemdoc, see http://jemdoc.jaboc.net/">
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<link rel="stylesheet" href="Non-negative%20Matrix%20Factorization%20under%20Heavy%20Noise_files/jemdoc.css" type="text/css">
</head>
<body>
<div id="layout-content">
<div id="toptitle">
<h1>Non-negative Matrix Factorization under Heavy Noise</h1>
</div>
<p><a href="http://drona.csa.iisc.ernet.in/%7Echiru/">Chiranjib Bhattacharyya</a>, <a href="http://research.microsoft.com/en-us/people/navingo/">Navin Goyal</a>, <a href="http://research.microsoft.com/en-us/people/kannan/">Ravindran Kannan</a>, Jagdeep Pani</a>.<br>
<a href="http://icml.cc/2016/">International Conference on Machine Learning (ICML), 2016</a>. <br><br></p>
<p><a href="http://mllab.csa.iisc.ernet.in/tsvdnmf/icml_main.pdf">Download Paper</a><br>
<a href="http://mllab.csa.iisc.ernet.in/tsvdnmf/icml_supp.pdf">Download Supplementary</a></p>
<h2>Abstract</h2>
<div class="infoblock">
<div class="blockcontent">
<p>The Noisy Non-negative Matrix factorization (NMF) is: given a data matrix <i>A</i> (<i>d</i> x <i>n</i>), find non-negative matrices <i>B</i>;<i>C</i> (<i>d</i> x <i>k</i>, <i>k</i> x <i>n</i> respy.) so that <i>A</i> = <i>BC</i> +<i>N</i>, where <i>N</i> is a noise matrix. Existing polynomial time algorithms with proven error guarantees require EACH column <i>N</i>⋅<i>j</i> to have <i>l1</i> norm much smaller than ||(<i>BC</i>)⋅<i>j</i>||1, which could be very restrictive. In important applications of NMF such as Topic Modeling as well as theoretical noise models (e.g. Gaussian with high sigma), almost EVERY column of <i>N</i>⋅<i>j</i> violates this condition. We introduce the heavy noise model which only requires the average noise over large subsets of columns to be small. We initiate a study of Noisy NMF under the heavy noise model. We show that our noise model subsumes noise models of theoretical and practical interest (for e.g. Gaussian noise of maximum possible sigma). We then devise an algorithm <i>TSVDNMF</i> which under certain assumptions on <i>B</i>,<i>C</i>, solves the problem under heavy noise. Our error guarantees match those of previous algorithms. Our running time of <i>O(k.(d+n)^2)</i> is substantially better than the <i>O(d.n^3)</i> for the previous best. Our assumption on <i>B</i> is weaker than the “Separability” assumption made by all previous results. We provide empirical justification for our assumptions on <i>C</i>. We also provide the first proof of identifiability (uniqueness of <i>B</i>) for noisy NMF which is not based on separability and does not use hard to check geometric conditions. Our algorithm outperforms earlier polynomial time algorithms both in time and error, particularly in the presence of high noise.</p>
</div></div>
<h2>How to run the Code</h2>
To be added.<br>
For any bugs, queries or suggestions, drop a mail to “pani dot jagdeep at gmail dot com”<p></p>
<h2>Contact</h2>
<p>You are welcome to contact the authors for any queries, thoughts or criticism. <br> Thanks ! (not factorial)
</p>

</body></html>
