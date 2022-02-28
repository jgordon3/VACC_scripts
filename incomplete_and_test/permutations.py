#/usr/local/python
#my_string = 'dakota'
#mapping = { 'd':'D', 'a':'A', 'k':'K', 'o':'0', 't':'T'}
#for k, v in mapping.iteritems():
#    my_string = my_string.replace(k, v)
#print my_string


#this is actually ruby
#stackoverflow.com/questions/22067184/all-possible-combinations-of-selected-character-substitution-in-a-string-in-ruby

string = '&dakota!'
subs = ['d'=>'D','a'=>'A','k'=>'K','o'=>'0','t'=>'T','a'=>'@','o'=>'O']

subs = subs.first.map(&:to_a)
1.upto(subs.length).each do |n|
    subs.combination(n).each do |a|
    p a.each_with_object(string.dup){|pair, s| s.gsub!(*pair)}
        end
end