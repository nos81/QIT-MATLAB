function U=two(g,targn,totn)
% TWO gives the unitary of a two qubit gate
%
% U=two(gate,num_targ_1,num_tot_q)
%
% Returns unitary gate for a two qubit gate applied to
% two sequential qubits and identity applied to the
% remaining qubits
%
% Counts of qubits begins at 1.
%
% Also see:
%     gate.controlled
%
% James Whitfield 2010

%TODO make bigger identity matrices to save time and code

%Warning
if size(g,1)~=4 | size(g,2)~=4
	warning('The gate is not a two qubit gate')
end


global qit;
id=qit.I;
eval_str='mkron(';
for k=1:(targn-1)
	eval_str=[eval_str,'id,'];
end
eval_str=[eval_str,'g'];
for k=(targn+2):totn
	eval_str=[eval_str,',id'];
end
eval_str=[eval_str,');'];

U=eval(eval_str);
