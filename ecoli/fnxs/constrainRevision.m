function feasible = constrainRevision(points,A,lb,ub,tol)
% Return true if point acomplish all the restricitons with <= tol
% tolerance.

lbFeas   = ~any(points-lb<-tol);
ubFeas   = ~any(points-ub>tol);
consFeas = ~any(abs(A * points) > tol);
feasible = lbFeas.*ubFeas.*consFeas;

end
