#!/bin/tcsh

set return_array = ()

set array = ()
foreach item ($4)
	set array = ($array $item)
end



set i=1
foreach item ($3)

	if ("$array[$i]" != $6) then

		if ($5 != 'plotOnly' && $5 != 'analyseTargetSelection' && $5 !~ plotXY_*) then
			if (! -e $1'myRun/plots/'$item'/') 	mkdir $1'myRun/plots/'$item'/'
			if (! -e $1'myRun/histos/'$item'/') 	mkdir $1'myRun/histos/'$item'/'
		endif

		touch $2'plot_'$item'_config.txt'
		rm -f $2'plot_'$item'_config.txt'

		set return_array=($return_array $2'plot_'$item'_config.txt')
	endif
	@ i++
end

echo $return_array
