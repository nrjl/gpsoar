function hyperin = positivify(hyperin, n_hyper, order_h)

n_nonmeanhyper = sum(n_hyper(1:end-1));
n_y = ((numel(hyperin)-n_nonmeanhyper)/(2*order_h));
indexes = zeros(1, n_y*order_h);
for i = 1:n_y
	indexes(((i-1)*order_h+1):(i*order_h)) = ((n_nonmeanhyper+1+(2*i-1)*order_h):(n_nonmeanhyper+i*order_h*2));
end

hyperin(indexes) = abs(hyperin(indexes));