function [G, Queue] = PopRotation(Queue)
%POPROTATION Pop a rotation from the queue. 

G = [ Queue(:,1) , [-conj(Queue(2,1)) ; conj(Queue(1,1)) ] ];
Queue = Queue(:,2:end);

end

