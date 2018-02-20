function topic = help_empty_topics(topic,k,~)
% This portion takes care of orphan topics i.e. it ensures that each topic has atleast one document
% save ('topic_temp_test.mat','topic','k','n');

for adjust_topic=1:k
    topic_count2(adjust_topic)=sum(topic==adjust_topic);
end
orphan_topics=find(topic_count2==0);
[~,max_top]=max(topic_count2);

doc_id2 = 1;
for all_count=1:length(orphan_topics)
    found1=0;
    while(~found1)
        if topic(doc_id2)==max_top
            found1=1;
            topic(doc_id2)=orphan_topics(all_count);
        end
        
        doc_id2=doc_id2+1;
    end
end

