function U=single(g,targn,totn)
% SINGLE gives the unitary for a single qubit gate
%
% U=single(gate,num_target_qubit,num_total_qubits)
%
% Returns unitary gate for a single qubit applied to target
% qubit and identity applied to the remaining qubits. 
%
% Counting of the qubits begins at 1. 
%
% James Whitfield 2010

%TODO make bigger identity matrices to save time and code

%Warnings
if targn>totn
	warning('the target qubit number is larger than the total number of qubits')
end
if size(g,1)~=2 | size(g,2)~=2
	warning('The gate is not a single qubit gate')
end


id=eye(2);

%create string for evaluation
eval_str='mkron(';
for k=1:(targn-1)
	eval_str=[eval_str,'id,'];
end
eval_str=[eval_str,'g'];
for k=(targn+1):totn
	eval_str=[eval_str,',id'];
end
eval_str=[eval_str,');'];

U=eval(eval_str);
