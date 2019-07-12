function y=Ids00_cMOSFET(QS,QD,Vth,Cox,Q0)
	ids1 = 2*Vth*(QS-QD)+(QS^2-QD^2)/Cox
	ids2 = Vth*Q0*ln((QD+Q0)/(QS+Q0))
	y=ids1+ids2
endfunction
