#! /usr/local/bin/perl5
#

for ($i = 0; $i < 1000; ++$i)
{
	print sin($i)+sin($i/500)+sin($i/100)+sin($i/25)+rand(0.5)*sin($i), "\n";
}



